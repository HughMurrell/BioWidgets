#!/usr/bin/env python3
"""
Annotate FASTA sequence identities with specimen collection dates.

Usage:
    python fasta_annotate.py <input.fasta> <samples.csv>

The sequence identity is the header token before the first space. It is split
on '_'. The PID is the first field with the prefix CAP replaced by 100
(CAP256 -> 100256). The visit code is the first four characters after that
underscore. Those two values are looked up in the CSV columns PID and
VisitCode. SpecimenCollectionDate values such as 16-May-19 are read as
years 2000-2099 and appended as YYYY-MM-DD:

    CAP256_2000AAAAACCG fs=44 ma=0.86
    -> CAP256_2000AAAAACCG|2008-03-15

Text after the first space (for example fs= and ma=) is dropped. Headers that
do not match this identity pattern are left unchanged. The output file is
written next to the input, with _annotated inserted before the extension
(input.fasta -> input_annotated.fasta).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path


REQUIRED_COLUMNS = ("PID", "VisitCode", "SpecimenCollectionDate")
MONTHS = {
    "jan": 1,
    "feb": 2,
    "mar": 3,
    "apr": 4,
    "may": 5,
    "jun": 6,
    "jul": 7,
    "aug": 8,
    "sep": 9,
    "oct": 10,
    "nov": 11,
    "dec": 12,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate FASTA sequence identities with specimen collection dates."
    )
    parser.add_argument("fasta", type=Path, help="Input FASTA file")
    parser.add_argument(
        "csv",
        type=Path,
        help="CSV with PID, VisitCode, and SpecimenCollectionDate columns",
    )
    return parser.parse_args()


def normalize_lookup_value(value: str) -> str:
    """Strip whitespace and a trailing .0 left by spreadsheet exports."""
    text = (value or "").strip()
    if len(text) > 2 and text.endswith(".0") and text[:-2].isdigit():
        return text[:-2]
    return text


def output_path_for(fasta_path: Path) -> Path:
    return fasta_path.with_name(f"{fasta_path.stem}_annotated{fasta_path.suffix}")


def format_collection_date(value: str) -> str | None:
    """Return a CSV collection date as YYYY-MM-DD, or None if it cannot be read."""
    text = value.strip()
    parsed = _parse_date(text)
    if parsed is not None:
        return parsed
    head = re.split(r"[T\s]", text, maxsplit=1)[0]
    if head != text:
        return _parse_date(head)
    return None


def _parse_date(text: str) -> str | None:
    iso = re.fullmatch(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", text)
    if iso:
        year, month, day = (int(part) for part in iso.groups())
        return _validated_date(year, month, day)

    # Day-first numeric dates, as used for South African specimen tables.
    numeric = re.fullmatch(r"(\d{1,2})[/-](\d{1,2})[/-](\d{4})", text)
    if numeric:
        day, month, year = (int(part) for part in numeric.groups())
        return _validated_date(year, month, day)

    # CSV dates look like 16-May-19. Two-digit years are 2000-2099.
    named = re.fullmatch(r"(\d{1,2})[- ]([A-Za-z]{3,9})[- ](\d{2}|\d{4})", text)
    if named:
        month = MONTHS.get(named.group(2)[:3].lower())
        if month is None:
            return None
        year = int(named.group(3))
        if year < 100:
            year += 2000
        return _validated_date(year, month, int(named.group(1)))
    return None


def _validated_date(year: int, month: int, day: int) -> str | None:
    try:
        return datetime(year, month, day).strftime("%Y-%m-%d")
    except ValueError:
        return None


def load_collection_dates(csv_path: Path) -> dict[tuple[str, str], str]:
    dates: dict[tuple[str, str], str] = {}
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SystemExit(f"Error: {csv_path} has no header row")

        field_map = {name.strip(): name for name in reader.fieldnames if name and name.strip()}
        missing = [column for column in REQUIRED_COLUMNS if column not in field_map]
        if missing:
            found = ", ".join(field_map) or "(none)"
            raise SystemExit(
                f"Error: {csv_path} is missing column(s): {', '.join(missing)}. Found: {found}"
            )

        for row_number, row in enumerate(reader, start=2):
            pid = normalize_lookup_value(row.get(field_map["PID"]) or "")
            visit = normalize_lookup_value(row.get(field_map["VisitCode"]) or "")
            date = (row.get(field_map["SpecimenCollectionDate"]) or "").strip()
            if not pid or not visit:
                print(
                    f"Warning: skipping CSV row {row_number} with a blank PID or VisitCode",
                    file=sys.stderr,
                )
                continue
            if not date:
                print(
                    f"Warning: CSV row {row_number} has no SpecimenCollectionDate "
                    f"for PID {pid} visit {visit}",
                    file=sys.stderr,
                )
                continue
            formatted = format_collection_date(date)
            if formatted is None:
                print(
                    f"Warning: CSV row {row_number} has an unrecognised "
                    f"SpecimenCollectionDate {date!r} for PID {pid} visit {visit}",
                    file=sys.stderr,
                )
                continue
            date = formatted

            key = (pid, visit)
            previous = dates.get(key)
            if previous is None:
                dates[key] = date
            elif previous != date:
                print(
                    f"Warning: PID {pid} visit {visit} has multiple collection dates "
                    f"({previous}, {date}); using {previous}",
                    file=sys.stderr,
                )
    return dates


def read_fasta(fasta_path: Path) -> list[tuple[str, list[str]]]:
    records: list[tuple[str, list[str]]] = []
    header: str | None = None
    sequence_lines: list[str] = []

    with fasta_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.rstrip("\n\r")
            if text.startswith(">"):
                if header is not None:
                    records.append((header, sequence_lines))
                header = text[1:]
                sequence_lines = []
            elif header is None:
                if text.strip():
                    raise SystemExit(
                        f"Error: {fasta_path} line {line_number} has sequence data before a header"
                    )
            else:
                sequence_lines.append(text)

    if header is not None:
        records.append((header, sequence_lines))
    if not records:
        raise SystemExit(f"Error: {fasta_path} contains no FASTA records")
    return records


def pid_and_visit(identity: str) -> tuple[str, str] | None:
    """Return (PID, visit code) for an identity such as CAP256_2000AAAAACCG."""
    pid_token, separator, remainder = identity.partition("_")
    if not separator or not pid_token.startswith("CAP") or len(remainder) < 4:
        return None
    return "100" + pid_token[3:], remainder[:4]


def annotate_header(header: str, dates: dict[tuple[str, str], str]) -> tuple[str, str | None]:
    """
    Return the rewritten header and an error message.

    An error message means the identity looked like a sample name but had no
    collection date. Other headers are returned unchanged and without an error.
    """
    identity = header.split(" ", 1)[0]
    keys = pid_and_visit(identity)
    if keys is None:
        return header, None

    date = dates.get(keys)
    if date is None:
        pid, visit = keys
        return identity, f"{identity} (PID {pid}, visit {visit})"
    return f"{identity}|{date}", None


def write_fasta(output_path: Path, records: list[tuple[str, list[str]]]) -> None:
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence_lines in records:
            handle.write(f">{header}\n")
            for sequence_line in sequence_lines:
                handle.write(f"{sequence_line}\n")


def main() -> int:
    args = parse_args()
    if not args.fasta.is_file():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        return 1
    if not args.csv.is_file():
        print(f"Error: CSV file not found: {args.csv}", file=sys.stderr)
        return 1

    dates = load_collection_dates(args.csv)
    records = read_fasta(args.fasta)
    annotated: list[tuple[str, list[str]]] = []
    unmatched: list[str] = []
    annotated_count = 0

    for header, sequence_lines in records:
        new_header, error = annotate_header(header, dates)
        if error is None and new_header != header:
            annotated_count += 1
        elif error is not None:
            unmatched.append(error)
        annotated.append((new_header, sequence_lines))

    output_path = output_path_for(args.fasta)
    write_fasta(output_path, annotated)

    print(f"Annotated {annotated_count} sequence(s)", flush=True)
    print(f"Wrote {output_path}", flush=True)
    if unmatched:
        print(
            f"Error: {len(unmatched)} sequence(s) had no SpecimenCollectionDate:",
            file=sys.stderr,
        )
        for identity in unmatched:
            print(f"  {identity}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
