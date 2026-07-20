from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, Iterator, Optional, List
from snippy_ng.exceptions import SnippyError


# SAM column indices (0-based)
QNAME, FLAG, RNAME, POS, CIGAR, SEQ = 0, 1, 2, 3, 5, 9

# FLAG bits
PAIRED = 0x1
MUNMAP = 0x8

_HAS_CLIP = re.compile(r"\d[SH]")
_START_CLIPS = re.compile(r"^(?:(\d+)H)?(?:(\d+)S)?")
_END_CLIPS_SUFFIX = re.compile(r"(?:(\d+)S)?(?:(\d+)H)?$")
_SO = re.compile(r"(?:^|\t)SO:([^\t\n\r]+)")


@dataclass(frozen=True)
class Clips:
    hl: int
    sl: int
    sr: int
    hr: int


class SamclipError(SnippyError):
    pass


def fai_to_dict(fai_lines: Iterable[str]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for line in fai_lines:
        if not line.strip():
            continue
        name, bp, *_ = line.rstrip("\n").split("\t")
        try:
            out[name] = int(bp)
        except ValueError:
            continue
    return out


def _parse_hd_so(header_lines: List[str]) -> str:
    for h in header_lines:
        if h.startswith("@HD"):
            m = _SO.search(h)
            if not m:
                raise SamclipError("Missing @HD SO tag; expected @HD ... SO:queryname")
            return m.group(1)
    raise SamclipError("Missing @HD header line; expected @HD ... SO:queryname")


def _split_header(sam_lines: Iterable[str]) -> tuple[List[str], Iterator[str]]:
    it = iter(sam_lines)
    header: List[str] = []
    try:
        for line in it:
            if line.startswith("@"):
                header.append(line)
            else:
                def records(first: str, rest: Iterator[str]) -> Iterator[str]:
                    yield first
                    yield from rest
                return header, records(line, it)
    except UnicodeDecodeError as e:
        raise SamclipError("Input is not valid UTF-8; cannot parse SAM header. Are you sure this is a SAM file?") from e
    return header, iter(())


def _end_clips(cigar: str) -> Clips:
    m_start = _START_CLIPS.match(cigar)
    m_end = _END_CLIPS_SUFFIX.search(cigar)
    return Clips(
        hl=int(m_start.group(1) or 0) if m_start else 0,
        sl=int(m_start.group(2) or 0) if m_start else 0,
        sr=int(m_end.group(1) or 0) if m_end else 0,
        hr=int(m_end.group(2) or 0) if m_end else 0,
    )


def _query_length(cigar: str) -> int:
    """Return the original query length, including hard-clipped bases."""
    return sum(
        int(length)
        for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if operation in {"M", "I", "S", "H", "=", "X"}
    )


def _bad(
    fields: List[str],
    contig_len: Dict[str, int],
    max_clip: Optional[int],
    max_clip_fraction: Optional[float],
) -> bool:
    cigar = fields[CIGAR]
    if not _HAS_CLIP.search(cigar):
        return False

    clips = _end_clips(cigar)

    try:
        start = int(fields[POS])
    except ValueError:
        raise SamclipError(f"Invalid POS: {fields[POS]!r}")

    seq = fields[SEQ]
    end = start + len(seq) - 1

    rname = fields[RNAME]
    clen = contig_len.get(rname)
    if clen is None:
        raise SamclipError(f"Reference contig '{rname}' not in contig lengths")

    L = clips.hl + clips.sl
    R = clips.hr + clips.sr
    if start <= 1 + L:
        L = 0
    if end >= clen - R:
        R = 0

    exceeds_length = max_clip is not None and ((L > max_clip) or (R > max_clip))
    clipped_bases = L + R
    query_length = _query_length(cigar)
    exceeds_fraction = (
        max_clip_fraction is not None
        and query_length > 0
        and (clipped_bases / query_length) > max_clip_fraction
    )
    return exceeds_length or exceeds_fraction


def samclip_filter_lines(
    sam_lines: Iterable[str],
    contig_lengths: Dict[str, int],
    *,
    max_clip: Optional[int] = 5,
    max_clip_fraction: Optional[float] = None,
    invert: bool = False,
    fix_mate: bool = True,
    on_debug: Optional[Callable[[str], None]] = None,
) -> Iterator[str]:
    """Filter SAM alignment lines based on terminal clipping.

    ``max_clip`` applies an absolute limit to either terminal end of an
    alignment after ignoring clipping that reaches the contig boundary. A read
    is filtered when either end exceeds the limit.

    ``max_clip_fraction`` applies a limit to the total terminal clipping as a
    proportion of the original query length, including hard-clipped bases. The
    numerator is the sum of terminal soft- and hard-clipping on both ends,
    again after boundary-aware clipping is removed.

    When both limits are set, an alignment failing either one is filtered.

    Optionally sets the MUNMAP flag on mates of filtered reads when appropriate.
    """
    if max_clip is not None and max_clip < 0:
        raise ValueError("max_clip must be >= 0")
    if max_clip_fraction is not None and not 0 <= max_clip_fraction <= 1:
        raise ValueError("max_clip_fraction must be between 0 and 1")

    header, records = _split_header(sam_lines)
    so = _parse_hd_so(header)
    if fix_mate and so != "queryname":
        raise SamclipError(f"Input must be queryname-sorted; header has SO:{so}. Use `samtools sort -n` to fix this.")

    for h in header:
        yield h

    if invert or not fix_mate:
        for line in records:
            # Optimize: partial split, avoid stripping newline
            fields = line.split("\t", 10)
            is_bad = _bad(fields, contig_lengths, max_clip, max_clip_fraction)

            if invert:
                if is_bad:
                    yield line
            else:
                if not is_bad:
                    yield line

        return

    group_q: Optional[str] = None
    group_lines: List[str] = []
    group_fields: List[List[str]] = []
    group_test: List[bool] = []

    def flush() -> Iterator[str]:
        nonlocal group_q, group_lines, group_fields, group_test
        if not group_lines:
            return

        any_bad = any(group_test)
        any_good = not all(group_test)
        do_fix_mate = (not invert) and fix_mate and any_bad and any_good

        for line, fields, is_bad in zip(group_lines, group_fields, group_test):
            if invert:
                if is_bad:
                    yield_line = line
                else:
                    continue
            else:
                if is_bad:
                    continue
                if do_fix_mate:
                    try:
                        flag = int(fields[FLAG])
                    except ValueError:
                        raise SamclipError(f"Invalid FLAG: {fields[FLAG]!r}")
                    if flag & PAIRED:
                        fields[FLAG] = str(flag | MUNMAP)
                        yield_line = "\t".join(fields)
                    else:
                        yield_line = line
                else:
                    yield_line = line

            if on_debug and do_fix_mate:
                on_debug(
                    f"QNAME={fields[QNAME]} keep={not is_bad} invert={invert} "
                    f"mate_fix={do_fix_mate}"
                )
            yield yield_line

        group_q, group_lines, group_fields, group_test = None, [], [], []
    
    for line in records:
        # Optimize: partial split, avoid stripping newline
        fields = line.split("\t", 10)
        q = fields[QNAME]

        if group_q is None:
            group_q = q
        if q != group_q:
            yield from flush()
            group_q = q

        group_lines.append(line)
        group_fields.append(fields)
        group_test.append(_bad(fields, contig_lengths, max_clip, max_clip_fraction))

    yield from flush()
