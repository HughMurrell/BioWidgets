#!/usr/bin/env python3
"""
Batch sequence-logo production (aliViz-compatible).

For each FASTA alignment in an input directory:
  1. Assume sequence 1 = reference, sequence 2 = founder (bottom / escape-from row)
  2. Group remaining sequences by name delimiter/field (defaults: '_' field 3)
  3. Map the named epitope (from CSV) onto reference coordinates → alignment columns
  4. Build a stacked escape logo per group (collapse-weighted when --collapse-field is set)
  5. Write escape_from_<epitope>_<alignment>.{png|svg|pdf}

Usage:
  python logo_production.py ALIGN_DIR EPITOPE_CSV EPITOPE OUT_DIR
  python logo_production.py ALIGN_DIR EPITOPE_CSV EPITOPE OUT_DIR --format svg
  python logo_production.py ALIGN_DIR EPITOPE_CSV EPITOPE OUT_DIR --palette ALIGNMENT
  python logo_production.py ALIGN_DIR EPITOPE_CSV EPITOPE OUT_DIR --delimiter _ --field 3
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
import tempfile
import xml.sax.saxutils as xmlesc
from collections import defaultdict
from pathlib import Path
from typing import Optional

# --- Constants matching aliViz.html ---
DEFAULT_DELIMITER = "_"
DEFAULT_FIELD = 3  # 1-based, same as aliViz Group dialog
DEFAULT_FORMAT = "png"
DEFAULT_PALETTE = "BIOCHEMICAL"
DEFAULT_POWER = 1.0
DEFAULT_CHAR_WIDTH = 15

FASTA_EXTENSIONS = {".fasta", ".fa", ".fas", ".fna", ".faa"}

# Logo UI labels → CONFIG.colorPalettes keys in aliViz
PALETTE_KEYS = {
    "ALIGNMENT": "Alignment (IUPAC)",
    "BIOCHEMICAL": "Biochemical",
}

COLOR_PALETTES: dict[str, dict[str, str]] = {
    "Alignment (IUPAC)": {
        "A": "#80a0f0", "R": "#f01505", "N": "#00ff00", "D": "#c048c0",
        "C": "#f08080", "Q": "#00ff00", "E": "#c048c0", "G": "#f09048",
        "H": "#15a4a4", "I": "#80a0f0", "L": "#80a0f0", "K": "#f01505",
        "M": "#80a0f0", "F": "#80a0f0", "P": "#ffff00", "S": "#00ff00",
        "T": "#00ff00", "W": "#80a0f0", "Y": "#15a4a4", "V": "#80a0f0",
        "*": "#999999", "-": "#ffffff", "default": "#ffffff",
    },
    "Biochemical": {
        "D": "#e74c3c", "E": "#e74c3c",
        "K": "#3498db", "R": "#3498db", "H": "#3498db",
        "S": "#2ecc71", "N": "#2ecc71", "T": "#2ecc71", "G": "#2ecc71", "Q": "#2ecc71",
        "F": "#f39c12", "A": "#f39c12", "P": "#f39c12", "I": "#f39c12", "L": "#f39c12",
        "V": "#f39c12", "M": "#f39c12", "W": "#f39c12", "Y": "#f39c12",
        "C": "#95a5a6", "-": "#bdc3c7", "*": "#999999", "default": "#ffffff",
    },
}

# Group / cluster palette (aliViz COLOR_PALETTE); group 0 → index 1
COLOR_PALETTE = [
    "#ff00ff", "#3b82f6", "#10b981", "#f59e0b", "#8b5cf6",
    "#06b6d4", "#84cc16", "#f97316", "#6366f1", "#14b8a6",
    "#64748b", "#fbbf24", "#a855f7", "#22c55e",
]

GAP_LOGO_FILL = "#6b7280"  # nabla / gap glyph
GAP_GLYPH = "∇"  # internal frequency key
BOTTOM_LABEL_FILL = "#ff00ff"
SEPARATOR_COLOR = "#d0d0d0"
FONT_AA = (
    "Courier New, Courier, &quot;DejaVu Sans&quot;, "
    "&quot;Arial Unicode MS&quot;, &quot;Segoe UI Symbol&quot;, sans-serif, monospace"
)


# ---------------------------------------------------------------------------
# FASTA / epitopes / grouping
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (sequence name, ungapped-name sequence string)."""
    records: list[tuple[str, str]] = []
    header: Optional[str] = None
    seq_parts: list[str] = []
    text = path.read_text(encoding="utf-8", errors="replace")
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts).upper()))
            header = line[1:].strip().split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def iter_alignment_files(directory: Path) -> list[Path]:
    files = [
        p for p in directory.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS
    ]
    return sorted(files, key=lambda p: p.name.lower())


def parse_epitopes_csv(path: Path) -> dict[str, list[dict[str, int]]]:
    """Parse aliViz epitope CSV: name,start:end|pos,... → {name: [{start,end}, ...]}."""
    epitopes: dict[str, list[dict[str, int]]] = {}
    text = path.read_text(encoding="utf-8", errors="replace")
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 2:
            continue
        name = parts[0]
        regions: list[dict[str, int]] = []
        for coord_str in parts[1:]:
            if not coord_str:
                continue
            if ":" in coord_str:
                a, b = coord_str.split(":", 1)
                try:
                    start, end = int(a.strip()), int(b.strip())
                except ValueError:
                    continue
                regions.append({"start": start, "end": end})
            else:
                try:
                    coord = int(coord_str)
                except ValueError:
                    continue
                regions.append({"start": coord, "end": coord})
        if regions:
            epitopes[name] = regions
    return epitopes


def build_ref_map(ref_seq: str) -> list[int]:
    """1-based running residue index on reference; gaps keep previous count (aliViz)."""
    ref_map: list[int] = []
    count = 0
    for ch in ref_seq:
        if ch != "-":
            count += 1
        ref_map.append(count)
    return ref_map


def columns_for_ref_region(
    ref_map: list[int], start: int, end: int
) -> list[int]:
    """Alignment columns for inclusive reference region [start, end] (aliViz)."""
    if not ref_map:
        return []
    first_col: Optional[int] = None
    for c, coord in enumerate(ref_map):
        if coord >= start:
            first_col = c
            break
    if first_col is None:
        return []
    last_col = first_col
    for c in range(first_col, len(ref_map)):
        if ref_map[c] > end:
            break
        last_col = c
    return list(range(first_col, last_col + 1))


def group_sequences(
    names: list[str],
    delimiter: str,
    field_1based: int,
    skip_indices: set[int],
) -> tuple[dict[str, int], dict[int, str], dict[int, list[int]]]:
    """
    Group by name field (aliViz). Returns:
      name→group_id, group_id→label, group_id→list of sequence indices
    Group IDs are 0..k-1 ordered by sorted unique labels.
    """
    field_idx = field_1based - 1
    labels_by_index: dict[int, str] = {}
    for i, name in enumerate(names):
        if i in skip_indices:
            continue
        parts = name.split(delimiter)
        if field_idx < 0 or field_idx >= len(parts):
            continue
        label = parts[field_idx]
        if label == "":
            continue
        labels_by_index[i] = label

    unique = sorted(set(labels_by_index.values()))
    label_to_gid = {lab: gid for gid, lab in enumerate(unique)}
    name_to_gid: dict[str, int] = {}
    gid_to_label: dict[int, str] = {gid: lab for lab, gid in label_to_gid.items()}
    members: dict[int, list[int]] = defaultdict(list)
    for i, label in labels_by_index.items():
        gid = label_to_gid[label]
        name_to_gid[names[i]] = gid
        members[gid].append(i)
    return name_to_gid, gid_to_label, dict(members)


def collapse_weight(name: str, delimiter: str, field_1based: Optional[int]) -> int:
    if field_1based is None:
        return 1
    parts = name.split(delimiter)
    idx = field_1based - 1
    if idx < 0 or idx >= len(parts):
        return 1
    try:
        w = int(parts[idx])
    except ValueError:
        return 1
    return w if w >= 1 else 1


def sanitize_filename(s: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_-]", "_", s)


def estimate_text_width(text: str, font_size: float) -> float:
    """Approximate bold sans-serif advance width (no DOM measurement)."""
    return len(text) * font_size * 0.56


def fit_font_size(text: str, available_width: float, start: float = 18.0, minimum: float = 8.0) -> float:
    size = start
    while size > minimum:
        if estimate_text_width(text, size) <= available_width:
            return size
        size -= 0.5
    return minimum


def aa_fill(letter: str, palette: dict[str, str]) -> str:
    if letter == GAP_GLYPH:
        return GAP_LOGO_FILL
    return palette.get(letter, palette.get("default", "#000000"))


def deepen_hex(color: str, factor: float = 0.82) -> str:
    """Darken a #RRGGBB colour toward black so SVG/PNG output looks closer to aliViz canvas."""
    if not color or not color.startswith("#") or len(color) < 7:
        return color
    try:
        r = int(color[1:3], 16)
        g = int(color[3:5], 16)
        b = int(color[5:7], 16)
    except ValueError:
        return color
    r = max(0, min(255, int(round(r * factor))))
    g = max(0, min(255, int(round(g * factor))))
    b = max(0, min(255, int(round(b * factor))))
    return f"#{r:02x}{g:02x}{b:02x}"


def bold_text_attrs(fill: str, stroke_width: float = 0.85) -> str:
    """Same-colour stroke makes glyph strokes read heavier under rsvg/cairo."""
    return (
        f'fill="{fill}" stroke="{fill}" stroke-width="{stroke_width}" '
        f'paint-order="stroke fill" stroke-linejoin="round"'
    )


# Cap-height fraction of font-size for Courier-like capitals (ink height, not
# the 0.68 stacking advance used between stacked letters in aliViz).
COURIER_CAP_RATIO = 0.52


def nabla_polygon_svg(
    cx: float,
    y_base: float,
    scale_y: float,
    base_font_size: float,
    fill: str,
    stroke_width: float = 0.6,
) -> str:
    """
    Draw nabla as an outlined (unfilled) triangle.

    Height matches Courier capital ink after the same vertical scale used for AA
    letters (`unit_cap_height * scale_y`). Drawn in final coordinates (no
    scale transform) so the stroke stays even.
    """
    if scale_y <= 0 or base_font_size <= 0:
        return ""
    # Same vertical sizing as AA: cap ink at font-size, then × scale_y
    height = base_font_size * COURIER_CAP_RATIO * scale_y
    # AA letters use scale(1, sy) so width is not vertically coupled; keep letter-like width
    width = base_font_size * 0.62
    half = width / 2.0
    top_y = y_base - height
    points = (
        f"{cx - half:.2f},{top_y:.2f} "
        f"{cx + half:.2f},{top_y:.2f} "
        f"{cx:.2f},{y_base:.2f}"
    )
    # Outline only — not a solid wedge; stroke scales gently with letter size
    sw = max(1.25, 0.55 * scale_y + 0.9)
    return (
        f'<polygon points="{points}" fill="none" stroke="{fill}" '
        f'stroke-width="{sw:.2f}" stroke-linejoin="round" stroke-linecap="round"/>'
    )


# ---------------------------------------------------------------------------
# Logo SVG (aliViz generateSequenceLogo)
# ---------------------------------------------------------------------------

def build_logo_svg(
    records: list[tuple[str, str]],
    file_name: str,
    epitope_name: str,
    regions: list[dict[str, int]],
    delimiter: str,
    field: int,
    palette_name: str,
    power: float = DEFAULT_POWER,
    char_width: int = DEFAULT_CHAR_WIDTH,
    collapse_field: Optional[int] = None,
) -> str:
    if len(records) < 2:
        raise ValueError("Alignment needs at least reference + founder (2 sequences)")

    ref_name, ref_seq = records[0]
    founder_name, founder_seq = records[1]
    bottom_index = 1
    bottom_seq = founder_seq
    bottom_label = f"Founder ({founder_name})"

    palette_key = PALETTE_KEYS.get(palette_name.upper(), "Biochemical")
    palette = COLOR_PALETTES[palette_key]

    # Pad sequences to equal length
    max_len = max(len(s) for _, s in records)
    seqs = [(n, s + ("-" * (max_len - len(s)))) for n, s in records]
    ref_seq = seqs[0][1]
    bottom_seq = seqs[bottom_index][1]
    names = [n for n, _ in seqs]

    ref_map = build_ref_map(ref_seq)

    # Exclude reference (index 0) and founder (index 1) from grouping.
    # Founder is the escape-from bottom row only, never a group logo.
    name_to_gid, gid_to_label, members = group_sequences(
        names, delimiter, field, skip_indices={0, bottom_index}
    )
    if not members:
        raise ValueError(
            f"No sequences could be grouped with delimiter={delimiter!r} field={field}"
        )

    group_ids = sorted(members.keys())

    region_column_ranges: list[dict] = []
    all_columns: set[int] = set()
    for region in regions:
        cols = columns_for_ref_region(ref_map, region["start"], region["end"])
        if not cols:
            continue
        region_column_ranges.append(
            {
                "start": min(cols),
                "end": max(cols),
                "region": region,
                "columnCount": len(cols),
            }
        )
        all_columns.update(cols)

    if not all_columns:
        raise ValueError(
            f"Epitope {epitope_name!r} maps to no alignment columns on this reference"
        )

    sorted_columns = sorted(all_columns)

    def weight_for(idx: int) -> int:
        return collapse_weight(names[idx], delimiter, collapse_field)

    group_frequencies: list[dict] = []
    for gid in group_ids:
        group_seqs_idx = list(members[gid])
        column_frequencies: dict[int, dict[str, float]] = {}

        for col in all_columns:
            bottom_aa = None
            bottom_is_gap = False
            if col < len(bottom_seq):
                aa = bottom_seq[col]
                if aa == "-":
                    bottom_is_gap = True
                else:
                    bottom_aa = aa

            freq_map: dict[str, float] = defaultdict(float)
            for i in group_seqs_idx:
                seq = seqs[i][1]
                if col >= len(seq):
                    continue
                w = weight_for(i)
                aa = seq[col]
                if aa == "-":
                        if not bottom_is_gap:
                            freq_map[GAP_GLYPH] += w
                else:
                    if bottom_is_gap or aa != bottom_aa:
                        freq_map[aa] += w
            if freq_map:
                column_frequencies[col] = dict(freq_map)

        weighted_count = sum(weight_for(i) for i in group_seqs_idx)
        group_frequencies.append(
            {
                "groupId": gid,
                "groupValue": gid_to_label[gid],
                "columnFrequencies": column_frequencies,
                "count": weighted_count,
            }
        )

    char_w = max(10, min(20, int(char_width)))
    logo_cell_height = 18.0
    base_font_size = 16.0
    base_height = base_font_size
    cap_height_ratio = 0.68
    logo_height = logo_cell_height
    group_logo_height = logo_height * 3.0
    pad_h = 40.0
    pad_v = 40.0
    pixels_per_char = 7

    max_label_width = len(bottom_label) * pixels_per_char if bottom_label else 0
    for gd in group_frequencies:
        gl = f"{gd['groupValue']} ({gd['count']})"
        max_label_width = max(max_label_width, len(gl) * pixels_per_char)
    left_label_width = max(150.0, max_label_width + 20.0)
    region_gap = 3.0

    num_groups = len(group_frequencies)
    n_region_gaps = max(0, len(region_column_ranges) - 1)
    svg_width = (
        len(sorted_columns) * char_w
        + n_region_gaps * region_gap
        + pad_h * 2
        + left_label_width
    )

    power_pct = (power - 1.0) * 100.0
    abs_pct = f"{abs(power_pct):.0f}"
    if abs(power_pct) < 0.5:
        # Default power=1.0: omit the zero-bias clause from the title
        title_text = f"Logogram for {epitope_name}, showing escape from Founder"
    elif power_pct < 0:
        title_text = (
            f"Logogram for {epitope_name}, showing escape from Founder, "
            f"with {abs_pct}% frequency bias towards rare bases"
        )
    else:
        title_text = (
            f"Logogram for {epitope_name}, showing escape from Founder, "
            f"with {abs_pct}% frequency bias towards abundant bases"
        )
    alignment_text = f"Alignment: {file_name}"
    available_width = svg_width - (pad_h + 10) * 2
    font_size = fit_font_size(title_text, available_width, 20.0)
    font_size2 = fit_font_size(alignment_text, available_width, font_size)
    title_height = font_size + 5 + font_size2 + 15
    svg_height = title_height + num_groups * group_logo_height + logo_height + pad_v * 2

    parts: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_width:.2f}" '
        f'height="{svg_height:.2f}" style="background:white;border:1px solid #ccc">',
        f'<text x="{pad_h + 10}" y="{pad_v + font_size}" font-family="sans-serif" '
        f'font-size="{font_size}" font-weight="900" fill="#000" '
        f'stroke="#000" stroke-width="0.4" paint-order="stroke fill">'
        f"{xmlesc.escape(title_text)}</text>",
        f'<text x="{pad_h + 10}" y="{pad_v + font_size + 5 + font_size2}" '
        f'font-family="sans-serif" font-size="{font_size2}" font-weight="900" fill="#000" '
        f'stroke="#000" stroke-width="0.3" paint-order="stroke fill">'
        f"{xmlesc.escape(alignment_text)}</text>",
    ]

    def add_hline(y: float, stroke_width: float = 1.0) -> None:
        parts.append(
            f'<line x1="{pad_h + left_label_width}" y1="{y}" '
            f'x2="{svg_width - pad_h}" y2="{y}" '
            f'stroke="{SEPARATOR_COLOR}" stroke-width="{stroke_width}"/>'
        )

    def add_vline(x: float, y1: float, y2: float, stroke_width: float = 2.0) -> None:
        parts.append(
            f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2}" '
            f'stroke="{SEPARATOR_COLOR}" stroke-width="{stroke_width}"/>'
        )

    if group_frequencies:
        top_sep_y = svg_height - pad_v - logo_height - (num_groups * group_logo_height)
        add_hline(top_sep_y)

    for i in range(0, num_groups + 1):
        if i == 0:
            sep_y = svg_height - pad_v
        else:
            sep_y = svg_height - pad_v - logo_height - ((i - 1) * group_logo_height)
        add_hline(sep_y)

    # Vertical region boundaries
    x_offset = pad_h + left_label_width
    y_top = title_height + pad_v
    y_bot = svg_height - pad_v
    add_vline(x_offset, y_top, y_bot)
    for i in range(len(region_column_ranges) - 1):
        r = region_column_ranges[i]
        x_offset += r["columnCount"] * char_w
        add_vline(x_offset, y_top, y_bot)
        x_offset += region_gap
    if region_column_ranges:
        end_x = (
            pad_h
            + left_label_width
            + len(sorted_columns) * char_w
            + n_region_gaps * region_gap
        )
        add_vline(end_x, y_top, y_bot)

    # Bottom founder row
    bottom_logo_y = svg_height - pad_v
    x_offset = pad_h + left_label_width
    region_idx = 0
    for col in sorted_columns:
        if region_idx < len(region_column_ranges) - 1:
            cur = region_column_ranges[region_idx]
            nxt = region_column_ranges[region_idx + 1]
            if col > cur["end"] and col == nxt["start"]:
                x_offset += region_gap
                region_idx += 1

        aa = bottom_seq[col] if col < len(bottom_seq) else ""
        is_gap = aa == "-"
        is_dot = aa == "."
        bottom_aa = None if (is_gap or is_dot or not aa) else aa
        cx = x_offset + char_w / 2
        if is_gap or is_dot or bottom_aa:
            if is_gap:
                fill = deepen_hex(GAP_LOGO_FILL)
                parts.append(
                    nabla_polygon_svg(
                        cx,
                        bottom_logo_y,
                        scale_y=1.0,
                        base_font_size=base_font_size,
                        fill=fill,
                        stroke_width=0.9,
                    )
                )
            elif is_dot:
                fill, ch = deepen_hex(aa_fill("default", palette)), "."
                parts.append(
                    f'<text x="{cx}" y="{bottom_logo_y}" text-anchor="middle" '
                    f'dominant-baseline="alphabetic" font-family="{FONT_AA}" '
                    f'font-size="16" font-weight="900" {bold_text_attrs(fill, 0.9)}>'
                    f"{xmlesc.escape(ch)}</text>"
                )
            else:
                fill, ch = deepen_hex(aa_fill(bottom_aa, palette)), bottom_aa  # type: ignore[arg-type]
                parts.append(
                    f'<text x="{cx}" y="{bottom_logo_y}" text-anchor="middle" '
                    f'dominant-baseline="alphabetic" font-family="{FONT_AA}" '
                    f'font-size="16" font-weight="900" {bold_text_attrs(fill, 0.9)}>'
                    f"{xmlesc.escape(ch)}</text>"
                )

        # Region boundary coordinate labels
        for r in region_column_ranges:
            boundary = None
            if col == r["start"] and ref_map[col] == r["region"]["start"]:
                boundary = r["region"]["start"]
            elif col == r["end"] and ref_map[col] == r["region"]["end"]:
                boundary = r["region"]["end"]
            if boundary is not None:
                ly = svg_height - pad_v + 20
                coord_fill = deepen_hex("#444444", 0.95)
                parts.append(
                    f'<text x="{cx}" y="{ly}" text-anchor="middle" '
                    f'dominant-baseline="middle" font-family="sans-serif" '
                    f'font-size="10" font-weight="900" '
                    f'{bold_text_attrs(coord_fill, 0.35)} '
                    f'transform="rotate(90 {cx} {ly})">{boundary}</text>'
                )
                break

        x_offset += char_w

    # Group rows
    for row_idx, gd in enumerate(group_frequencies):
        y_base = svg_height - pad_v - logo_height - (row_idx * group_logo_height)
        x_offset = pad_h + left_label_width
        region_idx = 0
        total = gd["count"] if gd["count"] > 0 else 1

        for col in sorted_columns:
            if region_idx < len(region_column_ranges) - 1:
                cur = region_column_ranges[region_idx]
                nxt = region_column_ranges[region_idx + 1]
                if col > cur["end"] and col == nxt["start"]:
                    x_offset += region_gap
                    region_idx += 1

            col_freq = gd["columnFrequencies"].get(col)
            if col_freq:
                sorted_aas = sorted(col_freq.items(), key=lambda kv: -kv[1])
                raw = [c / total for _, c in sorted_aas]
                variant_sum = sum(raw)
                bottom_freq = 1.0 - variant_sum
                frequencies = list(raw)
                if abs(power - 1.0) > 1e-12:
                    all_f = raw + [bottom_freq]
                    transformed = [math.pow(max(f, 0.0), power) for f in all_f]
                    s = sum(transformed) or 1.0
                    frequencies = [t / s for t in transformed[:-1]]

                y_cur = y_base
                cx = x_offset + char_w / 2
                for i, (aa, _count) in enumerate(sorted_aas):
                    freq = frequencies[i]
                    char_height = base_height * freq * 3.0
                    if char_height <= 0:
                        continue
                    # Visual stacking advance (aliViz); ink height of AAs is shorter
                    visual_h = cap_height_ratio * char_height
                    fill = deepen_hex(aa_fill(aa, palette))
                    scale_y = char_height / base_height
                    stroke_w = 0.9 / max(scale_y, 0.15)
                    if aa == GAP_GLYPH:
                        # Same scaleY transform as AA text; unit size = Courier cap ink
                        parts.append(
                            nabla_polygon_svg(
                                cx,
                                y_cur,
                                scale_y=scale_y,
                                base_font_size=base_font_size,
                                fill=fill,
                                stroke_width=stroke_w,
                            )
                        )
                    else:
                        parts.append(
                            f'<text x="{cx}" y="{y_cur}" text-anchor="middle" '
                            f'dominant-baseline="alphabetic" font-family="{FONT_AA}" '
                            f'font-size="{base_font_size}" font-weight="900" '
                            f'{bold_text_attrs(fill, stroke_w)} '
                            f'transform="translate({cx}, {y_cur}) scale(1, {scale_y}) '
                            f'translate({-cx}, {-y_cur})">{xmlesc.escape(aa)}</text>'
                        )
                    y_cur = y_cur - visual_h

            x_offset += char_w

        label_color = deepen_hex(COLOR_PALETTE[(gd["groupId"] + 1) % len(COLOR_PALETTE)], 0.88)
        label = f"{gd['groupValue']} ({gd['count']})"
        parts.append(
            f'<text x="{pad_h + 10}" y="{y_base}" text-anchor="start" '
            f'dominant-baseline="text-after-edge" font-family="sans-serif" '
            f'font-size="12" font-weight="900" {bold_text_attrs(label_color, 0.35)}>'
            f"{xmlesc.escape(label)}</text>"
        )

    parts.append(
        f'<text x="{pad_h + 10}" y="{bottom_logo_y}" text-anchor="start" '
        f'dominant-baseline="text-after-edge" font-family="sans-serif" '
        f'font-size="12" font-weight="900" '
        f'{bold_text_attrs(deepen_hex(BOTTOM_LABEL_FILL, 0.88), 0.35)}>'
        f"{xmlesc.escape(bottom_label)}</text>"
    )
    # Silence unused ref_name lint in some checkers
    _ = ref_name
    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Export (SVG / PNG / PDF)
# ---------------------------------------------------------------------------

def convert_with_rsvg(svg_path: Path, out_path: Path, fmt: str) -> bool:
    exe = shutil.which("rsvg-convert")
    if not exe:
        return False
    cmd = [exe, "-f", fmt, "-o", str(out_path), str(svg_path)]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return out_path.is_file()
    except (subprocess.CalledProcessError, OSError):
        return False


def convert_with_cairosvg(svg_text: str, out_path: Path, fmt: str) -> bool:
    try:
        import cairosvg  # type: ignore
    except ImportError:
        return False
    try:
        if fmt == "png":
            cairosvg.svg2png(bytestring=svg_text.encode("utf-8"), write_to=str(out_path))
        elif fmt == "pdf":
            cairosvg.svg2pdf(bytestring=svg_text.encode("utf-8"), write_to=str(out_path))
        else:
            return False
        return out_path.is_file()
    except Exception:
        return False


def write_logo_output(svg_text: str, out_path: Path, fmt: str) -> None:
    fmt = fmt.lower()
    if fmt == "svg":
        out_path.write_text(svg_text, encoding="utf-8")
        return

    with tempfile.TemporaryDirectory() as tmp:
        svg_path = Path(tmp) / "logo.svg"
        svg_path.write_text(svg_text, encoding="utf-8")
        if convert_with_rsvg(svg_path, out_path, fmt):
            return
        if convert_with_cairosvg(svg_text, out_path, fmt):
            return

    raise RuntimeError(
        f"Cannot write {fmt.upper()}: install librsvg (rsvg-convert) or the "
        f"cairosvg Python package. SVG generation succeeded; use --format svg "
        f"as a fallback."
    )


def process_alignment(
    path: Path,
    out_dir: Path,
    epitope_name: str,
    regions: list[dict[str, int]],
    fmt: str,
    palette: str,
    delimiter: str,
    field: int,
    power: float,
    collapse_field: Optional[int],
) -> Path:
    records = read_fasta(path)
    if len(records) < 2:
        raise ValueError("fewer than 2 sequences (need reference + founder)")

    svg_text = build_logo_svg(
        records=records,
        file_name=path.name,
        epitope_name=epitope_name,
        regions=regions,
        delimiter=delimiter,
        field=field,
        palette_name=palette,
        power=power,
        collapse_field=collapse_field,
    )

    stem = sanitize_filename(path.stem)
    ep = sanitize_filename(epitope_name)
    out_path = out_dir / f"escape_from_{ep}_{stem}.{fmt.lower()}"
    write_logo_output(svg_text, out_path, fmt)
    return out_path


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Batch aliViz-style sequence logos for every alignment in a directory."
    )
    ap.add_argument("align_dir", type=Path, help="Directory of FASTA alignments")
    ap.add_argument("epitope_csv", type=Path, help="Epitope CSV (aliViz format)")
    ap.add_argument("epitope", help="Epitope name to use from the CSV")
    ap.add_argument("output_dir", type=Path, help="Directory for logo outputs")
    ap.add_argument(
        "--format",
        choices=("png", "svg", "pdf"),
        default=DEFAULT_FORMAT,
        help=f"Output format (default {DEFAULT_FORMAT})",
    )
    ap.add_argument(
        "--palette",
        choices=("ALIGNMENT", "BIOCHEMICAL"),
        default=DEFAULT_PALETTE,
        help=f"AA colour palette (default {DEFAULT_PALETTE})",
    )
    ap.add_argument(
        "--delimiter",
        default=DEFAULT_DELIMITER,
        help=f"Group name delimiter (default {DEFAULT_DELIMITER!r})",
    )
    ap.add_argument(
        "--field",
        type=int,
        default=DEFAULT_FIELD,
        help=f"1-based name field for grouping (default {DEFAULT_FIELD})",
    )
    ap.add_argument(
        "--power",
        type=float,
        default=DEFAULT_POWER,
        help=f"Frequency power transform (default {DEFAULT_POWER}; aliViz slider)",
    )
    ap.add_argument(
        "--collapse-field",
        type=int,
        default=None,
        metavar="N",
        help="Optional 1-based name field holding integer collapse counts",
    )
    args = ap.parse_args(argv)

    align_dir: Path = args.align_dir
    epitope_csv: Path = args.epitope_csv
    out_dir: Path = args.output_dir

    if not align_dir.is_dir():
        print(f"Error: align_dir is not a directory: {align_dir}", file=sys.stderr)
        return 1
    if not epitope_csv.is_file():
        print(f"Error: epitope_csv not found: {epitope_csv}", file=sys.stderr)
        return 1

    epitopes = parse_epitopes_csv(epitope_csv)
    if args.epitope not in epitopes:
        available = ", ".join(sorted(epitopes)) or "(none)"
        print(
            f"Error: epitope {args.epitope!r} not in {epitope_csv.name}. "
            f"Available: {available}",
            file=sys.stderr,
        )
        return 1
    regions = epitopes[args.epitope]

    files = iter_alignment_files(align_dir)
    if not files:
        print(
            f"No FASTA files found in {align_dir} "
            f"(extensions: {sorted(FASTA_EXTENSIONS)})",
            file=sys.stderr,
        )
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)

    ok, fail = 0, 0
    for path in files:
        try:
            out = process_alignment(
                path,
                out_dir,
                epitope_name=args.epitope,
                regions=regions,
                fmt=args.format,
                palette=args.palette,
                delimiter=args.delimiter,
                field=args.field,
                power=args.power,
                collapse_field=args.collapse_field,
            )
            print(f"OK  {path.name} -> {out.name}")
            ok += 1
        except Exception as exc:
            print(f"FAIL {path.name}: {exc}", file=sys.stderr)
            fail += 1

    print(f"Done: {ok} ok, {fail} failed → {out_dir}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
