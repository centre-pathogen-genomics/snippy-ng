from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
import shutil


def _parse_info(info: str) -> dict[str, str]:
    if not info or info == ".":
        return {}
    parsed: dict[str, str] = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = ""
    return parsed


def _format_info(info: dict[str, str]) -> str:
    if not info:
        return "."
    return ";".join(key if value == "" else f"{key}={value}" for key, value in info.items())


def _is_snp(ref: str, alt: str) -> bool:
    alts = alt.split(",")
    return len(ref) == 1 and len(alts) == 1 and len(alts[0]) == 1 and not alts[0].startswith("<")


def _is_indel(ref: str, alt: str) -> bool:
    return any(not allele.startswith("<") and len(ref) != len(allele) for allele in alt.split(","))


def _nearest_distance(position: int, positions: list[int]) -> int | None:
    if not positions:
        return None
    index = bisect_right(positions, position)
    distances: list[int] = []
    if index < len(positions):
        distances.append(abs(positions[index] - position))
    if index > 0:
        distances.append(abs(position - positions[index - 1]))
    return min(distances) if distances else None


def _count_positions_in_window(position: int, positions: list[int], window: int) -> int:
    left = bisect_right(positions, position - window - 1)
    right = bisect_right(positions, position + window)
    return right - left


def _edge_distance(info: dict[str, str]) -> int | None:
    for key in ("CONTEXT_EDGE_DIST", "EDGE_DIST", "ASM_EDGE_DIST", "MUMMER_EDGE_DIST"):
        value = info.get(key)
        if value is None:
            continue
        try:
            return int(float(value))
        except ValueError:
            return None
    return None


def _append_filter(filter_value: str, label: str) -> str:
    if filter_value in {"", ".", "PASS"}:
        return label
    labels = filter_value.split(";")
    if label not in labels:
        labels.append(label)
    return ";".join(labels)


def _append_info_reason(info: dict[str, str], key: str, reason: str) -> None:
    current = info.get(key)
    if not current:
        info[key] = reason
        return
    reasons = current.split(",")
    if reason not in reasons:
        reasons.append(reason)
    info[key] = ",".join(reasons)


def filter_variant_context_vcf(
    input_vcf: Path,
    output_vcf: Path,
    max_local_snps: int = 0,
    local_snp_window: int = 0,
    min_snp_distance_to_indel: int = 0,
    min_snp_distance_to_breakpoint: int = 0,
) -> None:
    """Mark SNPs that fail enabled local-context rules as ``LowQual``."""
    local_snp_filter_enabled = max_local_snps > 0 and local_snp_window > 0
    indel_distance_filter_enabled = min_snp_distance_to_indel > 0
    breakpoint_distance_filter_enabled = min_snp_distance_to_breakpoint > 0

    if not (
        local_snp_filter_enabled
        or indel_distance_filter_enabled
        or breakpoint_distance_filter_enabled
    ):
        with open(input_vcf, "r", encoding="utf-8") as source, open(output_vcf, "w", encoding="utf-8") as destination:
            shutil.copyfileobj(source, destination)
        return

    headers: list[str] = []
    records: list[list[str]] = []
    snp_positions_by_ref: dict[str, list[int]] = defaultdict(list)
    indel_positions_by_ref: dict[str, list[int]] = defaultdict(list)

    with open(input_vcf, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                headers.append(raw_line)
                continue
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 8:
                records.append(fields)
                continue
            chrom, pos_text, _id, ref, alt = fields[:5]
            pos = int(pos_text)
            records.append(fields)
            if _is_snp(ref, alt):
                snp_positions_by_ref[chrom].append(pos)
            elif _is_indel(ref, alt):
                indel_positions_by_ref[chrom].append(pos)

    for positions in snp_positions_by_ref.values():
        positions.sort()
    for positions in indel_positions_by_ref.values():
        positions.sort()

    header_text = "".join(headers)
    context_header_lines = []

    context_header_lines.extend([
        ("FILTER", "LowQual", '##FILTER=<ID=LowQual,Description="Low-quality variant call based on local context">\n'),
        ("INFO", "CONTEXT_LOWQUAL_REASON", '##INFO=<ID=CONTEXT_LOWQUAL_REASON,Number=.,Type=String,Description="Context reasons this variant was marked LowQual">\n'),
    ])
    if indel_distance_filter_enabled:
        context_header_lines.append(
            ("INFO", "CONTEXT_INDEL_DIST", '##INFO=<ID=CONTEXT_INDEL_DIST,Number=1,Type=Integer,Description="Distance in reference bases from this SNP to the nearest indel candidate">\n')
        )
    if local_snp_filter_enabled:
        context_header_lines.append(
            ("INFO", "CONTEXT_LOCAL_SNPS", '##INFO=<ID=CONTEXT_LOCAL_SNPS,Number=1,Type=Integer,Description="Number of SNP candidates within the configured local SNP window">\n')
        )

    with open(output_vcf, "w", encoding="utf-8") as handle:
        for header in headers:
            if header.startswith("#CHROM"):
                for header_type, identifier, definition in context_header_lines:
                    if f"##{header_type}=<ID={identifier}," not in header_text:
                        handle.write(definition)
            handle.write(header)

        for fields in records:
            if len(fields) < 8:
                handle.write("\t".join(fields) + "\n")
                continue
            chrom, pos_text, _id, ref, alt = fields[:5]
            if not _is_snp(ref, alt):
                handle.write("\t".join(fields) + "\n")
                continue

            pos = int(pos_text)
            info = _parse_info(fields[7])
            lowqual_reasons: list[str] = []
            if breakpoint_distance_filter_enabled:
                edge_distance = _edge_distance(info)
                if edge_distance is not None and edge_distance < min_snp_distance_to_breakpoint:
                    lowqual_reasons.append("NEAR_ALIGNMENT_EDGE")
            if indel_distance_filter_enabled:
                indel_distance = _nearest_distance(pos, indel_positions_by_ref.get(chrom, []))
                if indel_distance is not None:
                    info["CONTEXT_INDEL_DIST"] = str(indel_distance)
                    if indel_distance < min_snp_distance_to_indel:
                        lowqual_reasons.append("NEAR_INDEL")
            if local_snp_filter_enabled:
                local_snps = _count_positions_in_window(pos, snp_positions_by_ref.get(chrom, []), local_snp_window)
                info["CONTEXT_LOCAL_SNPS"] = str(local_snps)
                if local_snps > max_local_snps:
                    lowqual_reasons.append("LOCAL_SNP_CLUSTER")
            if lowqual_reasons:
                for reason in lowqual_reasons:
                    _append_info_reason(info, "CONTEXT_LOWQUAL_REASON", reason)
                fields[6] = _append_filter(fields[6], "LowQual")
            fields[7] = _format_info(info)
            handle.write("\t".join(fields) + "\n")
