from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

from snippy_ng.exceptions import SnippyError


QNAME, FLAG, RNAME, POS, CIGAR, SEQ, QUAL = 0, 1, 2, 3, 5, 9, 10
UNMAP = 0x4

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_OPS = {"M", "D", "N", "=", "X"}
QUERY_OPS = {"M", "I", "S", "=", "X"}
QUERY_REF_OPS = {"M", "=", "X"}
TRIM_EDGE_OPS = {"D", "N", "P"}

Interval = Tuple[int, int]
IntervalsByContig = Dict[str, List[Interval]]


class SamcropError(SnippyError):
    pass


def parse_bed_lines(lines: Iterable[str]) -> IntervalsByContig:
    intervals: IntervalsByContig = defaultdict(list)
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            raise SamcropError(f"Invalid BED row: {line.rstrip()!r}")
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as error:
            raise SamcropError(f"Invalid BED coordinates: {line.rstrip()!r}") from error
        if start < 0 or end < start:
            raise SamcropError(f"Invalid BED interval: {line.rstrip()!r}")
        if end > start:
            intervals[fields[0]].append((start, end))
    return {contig: merge_intervals(values) for contig, values in intervals.items()}


def load_bed(path: Path) -> IntervalsByContig:
    with path.open("r", encoding="utf-8") as handle:
        return parse_bed_lines(handle)


def merge_intervals(intervals: Iterable[Interval]) -> List[Interval]:
    merged: List[Interval] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    if cigar == "*":
        return []
    parts = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not parts or "".join(f"{length}{op}" for length, op in parts) != cigar:
        raise SamcropError(f"Invalid CIGAR: {cigar!r}")
    return parts


def format_cigar(parts: Iterable[Tuple[int, str]]) -> str:
    compact: List[Tuple[int, str]] = []
    for length, op in parts:
        if length <= 0:
            continue
        if compact and compact[-1][1] == op:
            compact[-1] = (compact[-1][0] + length, op)
        else:
            compact.append((length, op))
    return "".join(f"{length}{op}" for length, op in compact)


def reference_length(cigar_parts: Iterable[Tuple[int, str]]) -> int:
    return sum(length for length, op in cigar_parts if op in REF_OPS)


def overlapping_intervals(read_start: int, read_end: int, intervals: List[Interval]) -> List[Interval]:
    overlaps: List[Interval] = []
    for start, end in intervals:
        if end <= read_start:
            continue
        if start >= read_end:
            break
        overlap_start = max(read_start, start)
        overlap_end = min(read_end, end)
        if overlap_start < overlap_end:
            overlaps.append((overlap_start, overlap_end))
    return overlaps


def _trim_edge_ops(parts: List[Tuple[int, str, Optional[int]]]) -> List[Tuple[int, str, Optional[int]]]:
    while parts and parts[0][1] in TRIM_EDGE_OPS:
        parts.pop(0)
    while parts and parts[-1][1] in TRIM_EDGE_OPS:
        parts.pop()
    return parts


def _first_reference_start(parts: List[Tuple[int, str, Optional[int]]], fallback: int) -> int:
    for _length, op, ref_start in parts:
        if op in REF_OPS and ref_start is not None:
            return ref_start
    return fallback


def _strip_stale_tags(fields: List[str]) -> List[str]:
    return fields[:11] + [
        field for field in fields[11:]
        if not (field.startswith("MD:Z:") or field.startswith("NM:i:") or field.startswith("SA:Z:"))
    ]


def _crop_sam_record_to_interval(
    fields: List[str],
    cigar_parts: List[Tuple[int, str]],
    read_start: int,
    keep: Interval,
) -> Optional[str]:
    fields = list(fields)
    keep_start, keep_end = keep

    seq = fields[SEQ]
    qual = fields[QUAL]
    has_seq = seq != "*"
    has_qual = qual != "*"
    ref_pos = read_start
    query_pos = 0
    leading_hard = 0
    trailing_hard = 0
    emitted_query = False
    retained_parts: List[Tuple[int, str, Optional[int]]] = []
    seq_parts: List[str] = []
    qual_parts: List[str] = []

    def hard_clip_query(length: int) -> None:
        nonlocal leading_hard, trailing_hard
        if length <= 0:
            return
        if emitted_query or retained_parts:
            trailing_hard += length
        else:
            leading_hard += length

    def append_part(length: int, op: str, ref_start: Optional[int]) -> None:
        if length <= 0:
            return
        retained_parts.append((length, op, ref_start))

    for length, op in cigar_parts:
        if op in QUERY_REF_OPS:
            op_ref_start = ref_pos
            op_ref_end = ref_pos + length
            keep_left = max(op_ref_start, keep_start)
            keep_right = min(op_ref_end, keep_end)
            left_query_clip = max(0, keep_left - op_ref_start)
            kept = max(0, keep_right - keep_left)
            right_query_clip = max(0, op_ref_end - keep_right)

            hard_clip_query(left_query_clip)
            if kept:
                query_start = query_pos + left_query_clip
                query_end = query_start + kept
                append_part(kept, op, keep_left)
                if has_seq:
                    seq_parts.append(seq[query_start:query_end])
                if has_qual:
                    qual_parts.append(qual[query_start:query_end])
                emitted_query = True
            hard_clip_query(right_query_clip)
            ref_pos += length
            query_pos += length
        elif op in {"D", "N"}:
            op_ref_start = ref_pos
            op_ref_end = ref_pos + length
            keep_left = max(op_ref_start, keep_start)
            keep_right = min(op_ref_end, keep_end)
            append_part(max(0, keep_right - keep_left), op, keep_left if keep_left < keep_right else None)
            ref_pos += length
        elif op == "I":
            if keep_start <= ref_pos < keep_end:
                append_part(length, op, ref_pos)
                if has_seq:
                    seq_parts.append(seq[query_pos:query_pos + length])
                if has_qual:
                    qual_parts.append(qual[query_pos:query_pos + length])
                emitted_query = True
            else:
                hard_clip_query(length)
            query_pos += length
        elif op == "S":
            hard_clip_query(length)
            query_pos += length
        elif op == "H":
            hard_clip_query(length)
        elif op == "P":
            if keep_start <= ref_pos < keep_end:
                append_part(length, op, None)
        else:
            raise SamcropError(f"Unsupported CIGAR operation: {op!r}")

    retained_parts = _trim_edge_ops(retained_parts)
    if not retained_parts or not any(op in QUERY_OPS for _length, op, _ref_start in retained_parts):
        return None

    new_cigar_parts: List[Tuple[int, str]] = []
    if leading_hard:
        new_cigar_parts.append((leading_hard, "H"))
    new_cigar_parts.extend((length, op) for length, op, _ref_start in retained_parts)
    if trailing_hard:
        new_cigar_parts.append((trailing_hard, "H"))

    fields[POS] = str(_first_reference_start(retained_parts, keep_start) + 1)
    fields[CIGAR] = format_cigar(new_cigar_parts)
    fields[SEQ] = "".join(seq_parts) if has_seq else "*"
    fields[QUAL] = "".join(qual_parts) if has_qual else "*"
    fields = _strip_stale_tags(fields)
    return "\t".join(fields) + "\n"


def crop_sam_records(line: str, intervals_by_contig: IntervalsByContig) -> List[str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        raise SamcropError(f"Invalid SAM record with fewer than 11 fields: {line.rstrip()!r}")

    try:
        flag = int(fields[FLAG])
        read_start = int(fields[POS]) - 1
    except ValueError as error:
        raise SamcropError(f"Invalid SAM FLAG/POS in record: {line.rstrip()!r}") from error

    if flag & UNMAP or fields[RNAME] == "*" or fields[CIGAR] == "*":
        return []

    cigar_parts = parse_cigar(fields[CIGAR])
    ref_len = reference_length(cigar_parts)
    if ref_len <= 0:
        return []

    read_end = read_start + ref_len
    overlaps = overlapping_intervals(read_start, read_end, intervals_by_contig.get(fields[RNAME], []))
    cropped_records = [
        cropped
        for keep in overlaps
        if (cropped := _crop_sam_record_to_interval(fields, cigar_parts, read_start, keep)) is not None
    ]
    return cropped_records


def crop_sam_record(line: str, intervals_by_contig: IntervalsByContig) -> Optional[str]:
    records = crop_sam_records(line, intervals_by_contig)
    return records[0] if records else None


def samcrop_filter_lines(
    sam_lines: Iterable[str],
    intervals_by_contig: IntervalsByContig,
) -> Iterator[str]:
    for line in sam_lines:
        if line.startswith("@"):
            yield line
            continue
        yield from crop_sam_records(line, intervals_by_contig)
