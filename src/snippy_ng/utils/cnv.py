from __future__ import annotations

import csv
import math
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, TextIO


class CNVError(ValueError):
    pass


@dataclass(frozen=True)
class CoverageRow:
    contig_id: str
    length: int
    read_depth: float


@dataclass(frozen=True)
class CNVRow:
    contig_id: str
    read_depth: float
    copy_number: int


@dataclass(frozen=True)
class FeatureRow:
    feature_id: str
    contig_id: str
    start: int
    end: int


@dataclass(frozen=True)
class FeatureCNVRow:
    feature_id: str
    contig_id: str
    start: int
    end: int
    read_depth: float
    copy_number: int


@dataclass(frozen=True)
class SingleCopyRegion:
    start: int
    end: int
    contig_id: str | None = None


def _round_half_up(value: float) -> int:
    return int(math.floor(value + 0.5))


def parse_samtools_coverage(lines: Iterable[str]) -> list[CoverageRow]:
    rows: list[CoverageRow] = []
    header: list[str] | None = None

    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue

        fields = line.split("\t")
        if fields[0].startswith("#"):
            header = [field.lstrip("#") for field in fields]
            continue

        if header:
            row = dict(zip(header, fields))
            contig_id = row.get("rname")
            endpos = row.get("endpos")
            read_depth = row.get("meandepth")
        else:
            if len(fields) < 7:
                raise CNVError(f"Invalid samtools coverage row: {line}")
            contig_id = fields[0]
            endpos = fields[2]
            read_depth = fields[6]

        if contig_id is None or endpos is None or read_depth is None:
            raise CNVError(f"Missing required samtools coverage fields in row: {line}")

        try:
            rows.append(
                CoverageRow(
                    contig_id=contig_id,
                    length=int(endpos),
                    read_depth=float(read_depth),
                )
            )
        except ValueError as exc:
            raise CNVError(f"Invalid numeric value in samtools coverage row: {line}") from exc

    if not rows:
        raise CNVError("No contig coverage rows found")

    return rows


def estimate_copy_numbers(coverage_rows: list[CoverageRow], baseline: float | None = None) -> list[CNVRow]:
    if not coverage_rows:
        raise CNVError("No contig coverage rows found")

    if baseline is None:
        baseline = baseline_depth(coverage_rows)
    elif baseline <= 0:
        raise CNVError("Baseline depth must be greater than zero")

    return [
        CNVRow(
            contig_id=row.contig_id,
            read_depth=row.read_depth,
            copy_number=_round_half_up(row.read_depth / baseline),
        )
        for row in coverage_rows
    ]


def baseline_depth(coverage_rows: list[CoverageRow]) -> float:
    if not coverage_rows:
        raise CNVError("No contig coverage rows found")

    baseline = max(coverage_rows, key=lambda row: row.length)
    if baseline.read_depth <= 0:
        raise CNVError(
            f"Largest contig '{baseline.contig_id}' has zero depth; cannot estimate copy numbers"
        )
    return baseline.read_depth


def largest_contig_id(coverage_rows: list[CoverageRow]) -> str:
    if not coverage_rows:
        raise CNVError("No contig coverage rows found")
    return max(coverage_rows, key=lambda row: row.length).contig_id


def write_cnv_table(rows: list[CNVRow], output: TextIO, include_header: bool = True) -> None:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    if include_header:
        writer.writerow(["contig_id", "read_depth", "copy_number"])
    for row in rows:
        writer.writerow([row.contig_id, f"{row.read_depth:g}", row.copy_number])


def parse_known_single_copy_region(region: str) -> SingleCopyRegion:
    region = region.strip()
    if not region:
        raise CNVError("Known single-copy region cannot be empty")

    named_match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
    if named_match:
        contig_id, start, end = named_match.groups()
        parsed = SingleCopyRegion(contig_id=contig_id, start=int(start), end=int(end))
    else:
        coord_match = re.fullmatch(r"(\d+),(\d+)", region)
        if not coord_match:
            raise CNVError(
                "Known single-copy region must be START,END or CONTIG:START-END"
            )
        start, end = coord_match.groups()
        parsed = SingleCopyRegion(start=int(start), end=int(end))

    if parsed.start < 1 or parsed.end < parsed.start:
        raise CNVError(f"Invalid known single-copy region: {region}")

    return parsed


def parse_gff_attributes(attributes: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in attributes.split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def parse_gff_features(gff: Path, feature_type: str = "CDS") -> list[FeatureRow]:
    features: list[FeatureRow] = []

    with gff.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                continue
            if fields[2] != feature_type:
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as exc:
                raise CNVError(f"Invalid GFF coordinates in row: {line}") from exc
            if start < 1 or end < start:
                raise CNVError(f"Invalid GFF interval in row: {line}")

            attributes = parse_gff_attributes(fields[8])
            feature_id = (
                attributes.get("ID")
                or attributes.get("Name")
                or attributes.get("gene")
                or attributes.get("locus_tag")
                or f"{fields[0]}:{start}-{end}"
            )
            features.append(
                FeatureRow(
                    feature_id=feature_id,
                    contig_id=fields[0],
                    start=start,
                    end=end,
                )
            )

    if not features:
        raise CNVError(f"No '{feature_type}' features found in {gff}")

    return features


def write_feature_bed(features: list[FeatureRow], output: TextIO) -> None:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    for feature in features:
        writer.writerow([feature.contig_id, feature.start - 1, feature.end, feature.feature_id])


def parse_samtools_depth(lines: Iterable[str]) -> dict[str, dict[int, int]]:
    depths: dict[str, dict[int, int]] = {}

    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            raise CNVError(f"Invalid samtools depth row: {line}")

        try:
            pos = int(fields[1])
            depth = int(fields[2])
        except ValueError as exc:
            raise CNVError(f"Invalid numeric value in samtools depth row: {line}") from exc

        depths.setdefault(fields[0], {})[pos] = depth

    return depths


def estimate_feature_copy_numbers(
    features: list[FeatureRow],
    depths: dict[str, dict[int, int]],
    baseline: float,
) -> list[FeatureCNVRow]:
    if baseline <= 0:
        raise CNVError("Baseline depth must be greater than zero")

    rows: list[FeatureCNVRow] = []
    for feature in features:
        contig_depths = depths.get(feature.contig_id, {})
        values = [
            contig_depths.get(pos, 0)
            for pos in range(feature.start, feature.end + 1)
        ]
        read_depth = float(median(values))
        rows.append(
            FeatureCNVRow(
                feature_id=feature.feature_id,
                contig_id=feature.contig_id,
                start=feature.start,
                end=feature.end,
                read_depth=read_depth,
                copy_number=_round_half_up(read_depth / baseline),
            )
        )

    return rows


def write_feature_cnv_table(
    rows: list[FeatureCNVRow],
    output: TextIO,
    include_header: bool = True,
) -> None:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    if include_header:
        writer.writerow(["feature_id", "contig_id", "start", "end", "read_depth", "copy_number"])
    for row in rows:
        writer.writerow(
            [
                row.feature_id,
                row.contig_id,
                row.start,
                row.end,
                f"{row.read_depth:g}",
                row.copy_number,
            ]
        )


def samtools_coverage(alignment: Path) -> str:
    command = ["samtools", "coverage", str(alignment)]

    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout


def samtools_depth(alignment: Path, bed: Path) -> str:
    command = ["samtools", "depth", "-aa", "-b", str(bed), str(alignment)]

    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout


def samtools_region_depth(alignment: Path, contig_id: str, start: int, end: int) -> str:
    command = ["samtools", "depth", "-aa", "-r", f"{contig_id}:{start}-{end}", str(alignment)]

    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout


def median_depth_for_region(
    alignment: Path,
    coverage_rows: list[CoverageRow],
    known_single_copy: str,
) -> float:
    region = parse_known_single_copy_region(known_single_copy)
    contig_id = region.contig_id or largest_contig_id(coverage_rows)
    depth_output = samtools_region_depth(
        alignment,
        contig_id=contig_id,
        start=region.start,
        end=region.end,
    )
    depths = parse_samtools_depth(depth_output.splitlines())
    contig_depths = depths.get(contig_id, {})
    values = [
        contig_depths.get(pos, 0)
        for pos in range(region.start, region.end + 1)
    ]
    baseline = float(median(values))
    if baseline <= 0:
        raise CNVError(
            f"Known single-copy region {contig_id}:{region.start}-{region.end} has zero median depth"
        )
    return baseline


def copy_number_variation(
    alignment: Path,
    gff: Path | None = None,
    feature_type: str = "gene",
    known_single_copy: str | None = None,
) -> list[CNVRow] | list[FeatureCNVRow]:
    coverage_output = samtools_coverage(alignment)
    coverage_rows = parse_samtools_coverage(coverage_output.splitlines())

    if known_single_copy:
        baseline = median_depth_for_region(
            alignment,
            coverage_rows,
            known_single_copy=known_single_copy,
        )
    else:
        baseline = baseline_depth(coverage_rows)

    if gff is None:
        return estimate_copy_numbers(coverage_rows, baseline=baseline)

    features = parse_gff_features(gff, feature_type=feature_type)
    with tempfile.NamedTemporaryFile("w", suffix=".bed") as bed_handle:
        write_feature_bed(features, bed_handle)
        bed_handle.flush()
        depth_output = samtools_depth(alignment, Path(bed_handle.name))

    depths = parse_samtools_depth(depth_output.splitlines())
    return estimate_feature_copy_numbers(features, depths, baseline)
