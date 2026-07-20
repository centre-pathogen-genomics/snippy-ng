"""Merge sorted VCFs and summarise their call-set overlap."""

from __future__ import annotations

import gzip
import html
import io
import subprocess
from collections import Counter
from pathlib import Path
from typing import Iterable, TextIO


class VcfMergeError(RuntimeError):
    """Raised when ``bcftools merge`` cannot complete."""


def merge_vcfs(input_vcfs: list[Path], output_vcf: Path | None = None, *, force_samples: bool = True) -> None:
    """Merge coordinate-sorted VCFs with bcftools without requiring indexes.

    ``force_samples`` keeps identically named samples as distinct columns, which
    is the usual case when comparing VCFs from alternate callers or runs.
    """
    if len(input_vcfs) < 2:
        raise VcfMergeError("At least two VCF files are required to merge")

    output_type = (
        "b"
        if output_vcf is not None and output_vcf.suffix == ".bcf"
        else "z"
        if output_vcf is not None and output_vcf.name.endswith(".vcf.gz")
        else "v"
    )
    command = [
        "bcftools",
        "merge",
        "--no-index",
        *( ["--force-samples"] if force_samples else [] ),
        "--output-type",
        output_type,
        *(["--output", str(output_vcf)] if output_vcf is not None else []),
        *map(str, input_vcfs),
    ]
    try:
        if output_vcf is not None:
            subprocess.run(command, check=True, capture_output=True, text=True)
        else:
            subprocess.run(command, check=True, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError as exc:
        raise VcfMergeError("bcftools was not found on PATH") from exc
    except subprocess.CalledProcessError as exc:
        message = (exc.stderr or exc.stdout or str(exc)).strip()
        raise VcfMergeError(message) from exc


def read_vcf_file_list(file_lists: Iterable[Path]) -> list[Path]:
    """Read one VCF path per non-empty, non-comment line from file lists."""
    paths: list[Path] = []
    for file_list in file_lists:
        for raw_line in file_list.read_text(encoding="utf-8").splitlines():
            value = raw_line.strip()
            if value and not value.startswith("#"):
                paths.append(Path(value).expanduser())
    return paths


def write_upset_svg(input_vcfs: list[Path], output_svg: Path, *, max_intersections: int = 40) -> None:
    """Write a compact UpSet plot showing per-file variant presence.

    Presence is based on canonical ``CHROM, POS, REF, ALT`` alleles, regardless
    of a record's FILTER value or genotype. Multiallelic records are split into
    one membership entry for each ALT allele.
    """
    svg = upset_svg_for_vcfs(input_vcfs, max_intersections=max_intersections)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_svg.write_text(svg, encoding="utf-8")


def upset_svg_for_vcfs(input_vcfs: list[Path], *, max_intersections: int = 40) -> str:
    names = _unique_names(input_vcfs)
    variant_sets = [_variant_keys(path) for path in input_vcfs]
    counts: Counter[tuple[bool, ...]] = Counter()
    for key in set().union(*variant_sets):
        counts[tuple(key in variants for variants in variant_sets)] += 1
    return _upset_svg(counts, names, max_intersections=max_intersections)


def upset_svg_from_merged_vcf_text(vcf_text: str, *, max_intersections: int = 40) -> str:
    if not vcf_text.strip():
        raise VcfMergeError("No VCF data received on stdin")
    return upset_svg_from_merged_vcf(io.StringIO(vcf_text), max_intersections=max_intersections)


def upset_svg_from_merged_vcf(handle: TextIO, *, max_intersections: int = 40) -> str:
    names: list[str] = []
    counts: Counter[tuple[bool, ...]] = Counter()

    for line in handle:
        if not line.strip():
            continue
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            fields = line.rstrip("\n").split("\t")
            names = [name.strip() for name in fields[9:] if name.strip()]
            continue
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        if not names:
            raise VcfMergeError("Merged VCF on stdin has no sample columns")

        format_fields = fields[8].split(":")
        try:
            gt_index = format_fields.index("GT")
        except ValueError:
            continue

        membership: list[bool] = []
        for sample_field in fields[9 : 9 + len(names)]:
            sample_values = sample_field.split(":")
            gt = sample_values[gt_index] if gt_index < len(sample_values) else "."
            membership.append(_gt_has_alt(gt))
        counts[tuple(membership)] += 1

    if not names:
        raise VcfMergeError("Merged VCF on stdin is missing a #CHROM header with sample names")

    return _upset_svg(counts, names, max_intersections=max_intersections)


def _gt_has_alt(gt: str) -> bool:
    if not gt or gt == ".":
        return False
    for allele in gt.replace("|", "/").split("/"):
        if allele.isdigit() and int(allele) > 0:
            return True
    return False


def _open_vcf(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(path, "r", encoding="utf-8")


def _variant_keys(path: Path) -> set[tuple[str, int, str, str]]:
    variants: set[tuple[str, int, str, str]] = set()
    with _open_vcf(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[4] in {".", "*"}:
                continue
            for alt in fields[4].split(","):
                if alt not in {".", "*"}:
                    variants.add((fields[0], int(fields[1]), fields[3], alt))
    return variants


def _unique_names(paths: list[Path]) -> list[str]:
    names: list[str] = []
    seen: Counter[str] = Counter()
    for path in paths:
        name = _vcf_sample_name(path) or path.name.removesuffix(".gz").removesuffix(".vcf").removesuffix(".bcf")
        seen[name] += 1
        names.append(name if seen[name] == 1 else f"{name} ({seen[name]})")
    return names


def _vcf_sample_name(path: Path) -> str | None:
    if path.suffix == ".bcf":
        return None
    try:
        with _open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) > 9 and fields[9].strip():
                        return fields[9].strip()
                    return None
    except OSError:
        return None
    return None


def _upset_svg(counts: Counter[tuple[bool, ...]], names: list[str], *, max_intersections: int) -> str:
    rows = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:max_intersections]
    if not rows:
        rows = [(tuple(False for _ in names), 0)]

    column_width, margin_left, margin_right, margin_top, bar_height = 30, 190, 36, 58, 250
    matrix_top, row_height = margin_top + bar_height + 36, 26
    width = margin_left + margin_right + len(rows) * column_width
    height = matrix_top + len(names) * row_height + 48
    maximum = max(count for _, count in rows) or 1
    set_sizes = [sum(count for membership, count in counts.items() if membership[index]) for index in range(len(names))]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif}.title{font-size:22px;font-weight:bold}.label{font-size:13px}.small{font-size:11px;fill:#475569}.count{font-size:11px;font-weight:bold}</style>',
        f'<rect width="{width}" height="{height}" fill="#fff"/>',
        '<text x="24" y="36" class="title">VCF call-set overlap</text>',
        f'<line x1="{margin_left}" y1="{margin_top + bar_height}" x2="{width - margin_right}" y2="{margin_top + bar_height}" stroke="#334155"/>',
    ]
    for tick in range(6):
        value = maximum * tick / 5
        y = margin_top + bar_height - value / maximum * bar_height
        lines.extend([
            f'<line x1="{margin_left}" y1="{y:.2f}" x2="{width - margin_right}" y2="{y:.2f}" stroke="#e2e8f0"/>',
            f'<text x="{margin_left - 8}" y="{y + 4:.2f}" text-anchor="end" class="small">{value:.0f}</text>',
        ])
    for index, (membership, count) in enumerate(rows):
        x = margin_left + index * column_width + column_width / 2
        height_value = count / maximum * bar_height
        y = margin_top + bar_height - height_value
        fill = "#2563eb" if all(membership) else "#dc2626"
        included = ", ".join(name for name, present in zip(names, membership) if present)
        lines.append(f'<rect x="{x - 7:.2f}" y="{y:.2f}" width="14" height="{height_value:.2f}" fill="{fill}"><title>{html.escape(included)}: {count}</title></rect>')
        lines.append(f'<text x="{x:.2f}" y="{y - 5:.2f}" text-anchor="middle" class="count">{count}</text>')
    for row_index, (name, size) in enumerate(zip(names, set_sizes)):
        y = matrix_top + row_index * row_height
        lines.append(f'<rect x="0" y="{y - 15}" width="{width}" height="{row_height}" fill="{"#f1f5f9" if row_index % 2 == 0 else "#fff"}"/>')
        lines.append(f'<text x="{margin_left - 50}" y="{y + 4}" text-anchor="end" class="label">{html.escape(name)}</text>')
        lines.append(f'<text x="{margin_left - 12}" y="{y + 4}" text-anchor="end" class="small">{size}</text>')
    for index, (membership, _count) in enumerate(rows):
        x = margin_left + index * column_width + column_width / 2
        active = [matrix_top + row_index * row_height for row_index, present in enumerate(membership) if present]
        if len(active) > 1:
            lines.append(f'<line x1="{x:.2f}" y1="{min(active):.2f}" x2="{x:.2f}" y2="{max(active):.2f}" stroke="#0f172a" stroke-width="2"/>')
        for row_index, present in enumerate(membership):
            y = matrix_top + row_index * row_height
            lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{5 if present else 3.5}" fill="{"#0f172a" if present else "#cbd5e1"}"/>')
    lines.append(f'<text x="{margin_left}" y="{height - 14}" class="small">Bars are canonical variants (CHROM, POS, REF, ALT); blue = shared by every VCF, red = other intersections. Left numbers are set sizes.</text>')
    return "\n".join(lines + ["</svg>"]) + "\n"
