#!/usr/bin/env python3
"""
Parse aliViz SVG filenames/title text and emit a CSV summary.

Extracted fields per SVG:
- CAP: CAP number from filename, e.g. CAP008 -> 8
- DPI: highest dpi+N value in filename
- CL: cluster count from filename token cl-N
- MaxDepth: value from SVG text like "Max depth: 0.00550"
- SEQS: value from SVG text like "s=79"

Usage:
    python EarlyMutationRateEstimate.py /path/to/svg_dir
    python EarlyMutationRateEstimate.py /path/to/svg_dir -o out.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path


CAP_RE = re.compile(r"CAP(\d+)", re.IGNORECASE)
DPI_RE = re.compile(r"dpi\+(\d+)", re.IGNORECASE)
CL_RE = re.compile(r"(?:^|[_-])cl-(\d+)(?:[_\.-]|$)", re.IGNORECASE)
MAX_DEPTH_RE = re.compile(r"Max depth:\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)
SEQS_RE = re.compile(r"\bs=(\d+)\b", re.IGNORECASE)


def parse_filename_fields(file_path: Path) -> dict[str, int | None]:
    name = file_path.name

    cap_match = CAP_RE.search(name)
    cap = int(cap_match.group(1)) if cap_match else None

    dpi_matches = [int(m) for m in DPI_RE.findall(name)]
    dpi = max(dpi_matches) if dpi_matches else None

    cl_match = CL_RE.search(name)
    cl = int(cl_match.group(1)) if cl_match else None

    return {"CAP": cap, "DPI": dpi, "CL": cl}


def parse_svg_text_fields(file_path: Path) -> dict[str, str | int | None]:
    text = file_path.read_text(encoding="utf-8", errors="replace")

    max_depth_match = MAX_DEPTH_RE.search(text)
    max_depth = max_depth_match.group(1) if max_depth_match else None

    seqs_match = SEQS_RE.search(text)
    seqs = int(seqs_match.group(1)) if seqs_match else None

    return {"MaxDepth": max_depth, "SEQS": seqs}


def build_record(file_path: Path) -> dict[str, str | int | None]:
    record: dict[str, str | int | None] = {
        "File": file_path.name,
    }
    record.update(parse_filename_fields(file_path))
    record.update(parse_svg_text_fields(file_path))
    return record


def collect_svg_files(directory: Path) -> list[Path]:
    return sorted(
        [p for p in directory.iterdir() if p.is_file() and p.suffix.lower() == ".svg"],
        key=lambda p: p.name.lower(),
    )


def write_csv(records: list[dict[str, str | int | None]], output_path: Path) -> None:
    fieldnames = ["File", "CAP", "DPI", "CL", "MaxDepth", "SEQS"]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def linear_regression(points: list[tuple[float, float]]) -> tuple[float, float] | None:
    if len(points) < 2:
        return None

    xs = [x for x, _ in points]
    ys = [y for _, y in points]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)

    sxx = sum((x - mean_x) ** 2 for x in xs)
    if sxx == 0:
        return None

    sxy = sum((x - mean_x) * (y - mean_y) for x, y in points)
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    return slope, intercept


def linear_regression_through_origin(points: list[tuple[float, float]]) -> tuple[float, float] | None:
    if len(points) < 2:
        return None

    denom = sum(x * x for x, _ in points)
    if denom == 0:
        return None

    slope = sum(x * y for x, y in points) / denom
    return slope, 0.0


def svg_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def format_float(value: float) -> str:
    if value == 0:
        return "0"
    abs_value = abs(value)
    if 1e-4 <= abs_value < 1e4:
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return f"{value:.6e}"


def write_scatter_svg(records: list[dict[str, str | int | None]], output_path: Path) -> int:
    points: list[dict[str, float | int | None]] = []
    for record in records:
        if record.get("CL") != 1:
            continue
        dpi = record.get("DPI")
        max_depth = record.get("MaxDepth")
        if dpi is None or max_depth is None:
            continue
        try:
            x = float(dpi)
            y = float(max_depth)
        except (TypeError, ValueError):
            continue
        if x > 175:
            continue
        points.append({"x": x, "y": y, "cap": record.get("CAP")})

    width = 900
    height = 640
    margin_left = 90
    margin_right = 40
    margin_top = 70
    margin_bottom = 80
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom

    title = "Early Mutation Rate Estimate"
    subtitle = f"CL=1 records: {len(points)}"

    if not points:
        svg = f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="{width}" height="{height}" fill="white"/>
  <text x="{width/2}" y="40" text-anchor="middle" font-family="sans-serif" font-size="24" font-weight="700">{svg_escape(title)}</text>
  <text x="{width/2}" y="80" text-anchor="middle" font-family="sans-serif" font-size="16" fill="#4b5563">No CL=1 records with both DPI and MaxDepth were found.</text>
</svg>
"""
        output_path.write_text(svg, encoding="utf-8")
        return 0

    xs = [float(p["x"]) for p in points]
    ys = [float(p["y"]) for p in points]
    x_min = 0.0
    x_max = max(xs)
    y_min = 0.0
    y_max = max(ys)

    if x_min == x_max:
        x_max += 1
    else:
        x_pad = 0.05 * (x_max - x_min)
        x_max += x_pad

    if y_min == y_max:
        y_max += max(0.0001, 0.05 * abs(y_max) if y_max else 0.0001)
    else:
        y_pad = 0.08 * (y_max - y_min)
        y_max += y_pad

    def x_to_px(x: float) -> float:
        return margin_left + ((x - x_min) / (x_max - x_min)) * plot_width

    def y_to_px(y: float) -> float:
        return margin_top + plot_height - ((y - y_min) / (y_max - y_min)) * plot_height

    xy_points = [(float(p["x"]), float(p["y"])) for p in points]
    regression = linear_regression_through_origin(xy_points)

    svg_parts: list[str] = [
        f'<?xml version="1.0" encoding="UTF-8"?>\n',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n',
        f'  <rect width="{width}" height="{height}" fill="white"/>\n',
        f'  <text x="{width/2}" y="34" text-anchor="middle" font-family="sans-serif" font-size="24" font-weight="700">{svg_escape(title)}</text>\n',
        f'  <text x="{width/2}" y="58" text-anchor="middle" font-family="sans-serif" font-size="14" fill="#4b5563">{svg_escape(subtitle)}</text>\n',
        f'  <rect x="{margin_left}" y="{margin_top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#d1d5db" stroke-width="1"/>\n',
    ]

    x_ticks = 6
    y_ticks = 6
    for i in range(x_ticks + 1):
        x_val = x_min + (i / x_ticks) * (x_max - x_min)
        x_px = x_to_px(x_val)
        svg_parts.append(
            f'  <line x1="{x_px:.2f}" y1="{margin_top}" x2="{x_px:.2f}" y2="{margin_top + plot_height}" stroke="#f3f4f6" stroke-width="1"/>\n'
        )
        svg_parts.append(
            f'  <text x="{x_px:.2f}" y="{margin_top + plot_height + 24}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#374151">{svg_escape(format_float(x_val))}</text>\n'
        )

    for i in range(y_ticks + 1):
        y_val = y_min + (i / y_ticks) * (y_max - y_min)
        y_px = y_to_px(y_val)
        svg_parts.append(
            f'  <line x1="{margin_left}" y1="{y_px:.2f}" x2="{margin_left + plot_width}" y2="{y_px:.2f}" stroke="#f3f4f6" stroke-width="1"/>\n'
        )
        svg_parts.append(
            f'  <text x="{margin_left - 12}" y="{y_px + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="12" fill="#374151">{svg_escape(format_float(y_val))}</text>\n'
        )

    svg_parts.extend(
        [
            f'  <line x1="{margin_left}" y1="{margin_top + plot_height}" x2="{margin_left + plot_width}" y2="{margin_top + plot_height}" stroke="#111827" stroke-width="1.5"/>\n',
            f'  <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_height}" stroke="#111827" stroke-width="1.5"/>\n',
            f'  <text x="{margin_left + plot_width / 2}" y="{height - 20}" text-anchor="middle" font-family="sans-serif" font-size="16">DPI</text>\n',
            f'  <text x="24" y="{margin_top + plot_height / 2}" text-anchor="middle" font-family="sans-serif" font-size="16" transform="rotate(-90 24 {margin_top + plot_height / 2})">MaxDepth</text>\n',
        ]
    )

    for point in points:
        x = float(point["x"])
        y = float(point["y"])
        svg_parts.append(
            f'  <circle cx="{x_to_px(x):.2f}" cy="{y_to_px(y):.2f}" r="4.5" fill="#2563eb" fill-opacity="0.85"/>\n'
        )

    if regression is not None:
        slope, intercept = regression
        x1 = x_min
        x2 = x_max
        y1 = slope * x1 + intercept
        y2 = slope * x2 + intercept
        svg_parts.append(
            f'  <line x1="{x_to_px(x1):.2f}" y1="{y_to_px(y1):.2f}" x2="{x_to_px(x2):.2f}" y2="{y_to_px(y2):.2f}" stroke="#dc2626" stroke-width="2.5"/>\n'
        )
        annotation = f"slope = {format_float(slope)}, intercept = 0"

        residuals = []
        for point in points:
            x = float(point["x"])
            y = float(point["y"])
            predicted = slope * x + intercept
            residuals.append(y - predicted)

        outlier_indices: set[int] = set()
        if residuals:
            mean_abs_residual = sum(abs(r) for r in residuals) / len(residuals)
            outlier_threshold = 2.0 * mean_abs_residual
            if outlier_threshold > 0:
                for idx, (point, residual) in enumerate(zip(points, residuals)):
                    if abs(residual) < outlier_threshold:
                        continue
                    outlier_indices.add(idx)
                    cap = point.get("cap")
                    if cap is None:
                        continue
                    x = float(point["x"])
                    y = float(point["y"])
                    label = f"CAP{int(cap):03d}"
                    svg_parts.append(
                        f'  <text x="{x_to_px(x) + 8:.2f}" y="{y_to_px(y) - 8:.2f}" text-anchor="start" font-family="sans-serif" font-size="12" fill="#111827">{svg_escape(label)}</text>\n'
                    )

        filtered_points = [
            (float(point["x"]), float(point["y"]))
            for idx, point in enumerate(points)
            if idx not in outlier_indices
        ]
        filtered_regression = linear_regression_through_origin(filtered_points)
        if filtered_regression is not None:
            filtered_slope, filtered_intercept = filtered_regression
            fy1 = filtered_slope * x1 + filtered_intercept
            fy2 = filtered_slope * x2 + filtered_intercept
            svg_parts.append(
                f'  <line x1="{x_to_px(x1):.2f}" y1="{y_to_px(fy1):.2f}" x2="{x_to_px(x2):.2f}" y2="{y_to_px(fy2):.2f}" stroke="#059669" stroke-width="2.5" stroke-dasharray="8 6"/>\n'
            )
            annotation += f" | no-outlier slope = {format_float(filtered_slope)}"
    else:
        annotation = "Regression unavailable (need at least 2 CL=1 points with varying DPI)."

    svg_parts.append(
        f'  <text x="{margin_left + 8}" y="{margin_top + 20}" text-anchor="start" font-family="sans-serif" font-size="13" fill="#111827">{svg_escape(annotation)}</text>\n'
    )
    svg_parts.append("</svg>\n")

    output_path.write_text("".join(svg_parts), encoding="utf-8")
    return len(points)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract CAP/DPI/CL/MaxDepth/SEQS from SVG files and write a CSV."
    )
    parser.add_argument(
        "directory",
        type=Path,
        help="Directory containing SVG files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output CSV path. Default: <directory>/EarlyMutationRateEstimate.csv",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    directory: Path = args.directory

    if not directory.exists():
        print(f"Error: directory does not exist: {directory}", file=sys.stderr)
        return 1
    if not directory.is_dir():
        print(f"Error: not a directory: {directory}", file=sys.stderr)
        return 1

    svg_files = collect_svg_files(directory)
    if not svg_files:
        print(f"Error: no SVG files found in {directory}", file=sys.stderr)
        return 1

    records = [build_record(svg_file) for svg_file in svg_files]
    output_path = args.output or (directory / "EarlyMutationRateEstimate.csv")
    write_csv(records, output_path)
    plot_output_path = directory / "EarlyMutationRateEstimate.svg"
    plotted_points = write_scatter_svg(records, plot_output_path)

    print(f"Wrote {len(records)} record(s) to {output_path}")
    print(f"Wrote CL=1 plot with {plotted_points} point(s) to {plot_output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
