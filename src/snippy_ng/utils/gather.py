from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union
from snippy_ng.exceptions import SnippyError
import gzip
import io
import re

class GatherSamplesError(SnippyError):
    pass


def _format_file_list(paths: Iterable[Path]) -> str:
    """Format a list of paths as bullet points."""
    return '\n'.join(f"- {str(p)}" for p in paths)


def _raise_multiple_files_error(sample_id: str, sample_type: str, paths: List[Path]) -> None:
    """Raise error for samples that should have exactly one file but have multiple."""
    raise GatherSamplesError(
        f"{sample_id}: {sample_type} sample has multiple files!\n\n"
        f"Found:\n{_format_file_list(paths)}"
    )


def _raise_duplicate_candidate_error(file_type: str, first: Path, second: Path) -> None:
    """Raise error when multiple candidates are found for the same file type."""
    raise GatherSamplesError(
        f"Multiple {file_type} candidates!\n\n"
        f"Found:\n- {str(first)}\n- {str(second)}"
    )

SeqKind = str  # "ILL" | "ONT" | "ASM" | ...

@dataclass(frozen=True)
class SeqFile:
    path: Path
    kind: SeqKind           # e.g. "ILL", "ONT", "ASM"
    sample_id: str          # inferred ID from filename


def _open_maybe_compressed(p: Path) -> io.TextIOBase:
    """
    Open text stream, handling gzip by extension.
    Extend here if you want xz/zstd later.
    """
    if p.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return p.open("r", encoding="utf-8", errors="replace")


def detect_seq_kind(path: Union[str, Path]) -> Optional[SeqKind]:
    """
    Peek first line:
      - FASTA: '>'
      - FASTQ: '@...'
         heuristic: ONT ids often look UUID-like or contain key=value tokens
         otherwise ILL
    Returns "ASM", "ONT", "ILL", or None.
    """
    p = Path(path)
    if not p.is_file():
        return None
    try:
        if p.stat().st_size <= 0:
            return None
    except OSError:
        return None

    try:
        with _open_maybe_compressed(p) as fh:
            first = fh.readline()
    except Exception:
        return None

    if not first:
        return None
    if first.startswith(">"):
        return "ASM"
    if first.startswith("@"):
        m = re.match(r"^@(\S+)", first)
        if not m:
            return None
        rid = m.group(1)

        # ONT heuristic (tune as you like)
        is_uuid = bool(re.match(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", rid, re.I))
        has_kv = bool(re.match(r"^\w+=\w+", rid))
        return "ONT" if (is_uuid or has_kv) else "ILL"

    return None


def guess_sample_id(filename: str, aggressive: bool = False) -> str:
    """
    Control-A1_S1_L001_R1_001.fastq.gz -> Control-A1_S1_L001
    If aggressive=True, also strip _S1 and _L001, leaving just Control-A1
    """
    name = filename
    # strip compression
    name = re.sub(r"\.(?:gz|xz|zstd)$", "", name, flags=re.I)
    # strip seq format
    name = re.sub(r"\.(?:fq|fastq|fa|fasta|fna)$", "", name, flags=re.I)

    # optional illumina lane / sample indices
    if aggressive:
        name = re.sub(r"_S\d{1,3}", "", name)
        name = re.sub(r"_L\d{3}", "", name)

    # strip common read-direction suffixes
    # examples: _R1, _R1_001, _1, _2
    name = re.sub(r"_(?:R?[12])(?:_\d+)?$", "", name)

    return name


def scan_sequence_files(
    inputs: Iterable[Union[str, Path]],
    max_depth: int = 4,
    aggressive_ids: bool = False,
    exclude_name_regex: Optional[str] = None,
    exclude_files: Optional[List[Path]] = None,
) -> List[SeqFile]:
    """
    Walk the given files/folders (recursively up to max_depth) and detect sequence files.
    exclude_name_regex is applied to the basename only (not full path).
    """
    exclude_re = re.compile(exclude_name_regex) if exclude_name_regex else None

    found: List[SeqFile] = []
    for inp in inputs:
        root = Path(inp).expanduser()
        if root.is_file():
            candidates = [root]
        elif root.is_dir():
            # manual depth control
            base_parts = len(root.resolve().parts)
            candidates = []
            for p in root.rglob("*"):
                p = p.resolve()
                try:
                    if not p.is_file():
                        continue
                    depth = len(p.parts) - base_parts
                    if depth > max_depth:
                        continue
                    candidates.append(p)
                except OSError:
                    continue
        else:
            continue

        for p in candidates:
            if exclude_files and p in exclude_files:
                continue
            if exclude_re and exclude_re.search(p.name):
                continue
            kind = detect_seq_kind(p)
            if not kind:
                continue
            sid = guess_sample_id(p.name, aggressive=aggressive_ids)
            found.append(SeqFile(path=p, kind=kind, sample_id=sid))

    return found


def _find_r1_r2(files: List[Path]) -> Tuple[Optional[Path], Optional[Path]]:
    r1 = None
    r2 = None
    for p in files:
        n = p.name
        if re.search(r"(?:_R?1(?:_\d+)?)(?:\.\w+)*$", n):
            if r1 is not None:
                _raise_duplicate_candidate_error("R1", r1, p)
            r1 = p
        elif re.search(r"(?:_R?2(?:_\d+)?)(?:\.\w+)*$", n):
            if r2 is not None:
                _raise_duplicate_candidate_error("R2", r2, p)
            r2 = p
    return r1, r2


def handle_ILL(sample_id: str, items: List[SeqFile]) -> Dict:
    # "short": left/right
    paths = sorted([x.path for x in items])
    r1, r2 = _find_r1_r2(paths)

    # Fallback if names don't encode R1/R2: take first two in sort order
    if r1 is None or r2 is None:
        if len(paths) == 2:
            r1, r2 = paths
        else:
            raise GatherSamplesError(
                f"{sample_id}: ILL sample needs 2 files (R1/R2)!\n\n"
                f"Found:\n{_format_file_list(paths)}"
            )

    return {
        "type": "short",
        "left": str(r1),
        "right": str(r2),
    }


def handle_ONT(sample_id: str, items: List[SeqFile]) -> Dict:
    paths = sorted([x.path for x in items])
    if not paths:
        raise GatherSamplesError(f"{sample_id}: ONT sample has no files")
    if len(paths) > 1:
        _raise_multiple_files_error(sample_id, "ONT", paths)
    entry = {
        "type": "long",
        "reads": str(paths[0]),
        "caller": "clair3",
        "clair3_model": "",
    }
    return entry


def handle_ASM(sample_id: str, items: List[SeqFile]) -> Dict:
    paths = sorted([x.path for x in items])
    if not paths:
        raise GatherSamplesError(f"{sample_id}: ASM sample has no files")
    if len(paths) > 1:
        _raise_multiple_files_error(sample_id, "ASM", paths)
    return {
        "type": "asm",
        "assembly": str(paths[0]),
    }


SampleHandler = Callable[[str, List[SeqFile], dict], Dict]

def build_samples_config(
    seqfiles: List[SeqFile],
) -> Dict:
    """
    Build a dict shaped like:
      {"samples": { "sample1": {...}, "sample2": {...} } }
    """
    handlers: Dict[SeqKind, SampleHandler] = {
        "ILL": handle_ILL,
        "ONT": handle_ONT,
        "ASM": handle_ASM,
    }
    # group by sample_id
    grouped: Dict[str, List[SeqFile]] = {}
    for sf in seqfiles:
        grouped.setdefault(sf.sample_id, []).append(sf)

    samples: Dict[str, Dict] = {}
    for sid, items in sorted(grouped.items()):
        kinds = {x.kind for x in items}
        if len(kinds) != 1:
            offending_files = '\n'.join([f"- {str(x.path)} ({x.kind})" for x in items])
            raise GatherSamplesError(
                f"Multiple kinds detected for sample_id '{sid}'.\n\n"
                f"Filename parsing must not be ambiguous. "
                f"Please rename files or consider adjusting aggressive_ids and exclude_name_regex settings.\n\n"
                f"Found:\n{offending_files}"
            )
        kind = next(iter(kinds))
        handler = handlers.get(kind)
        if handler is None:
            raise KeyError(
                f"{sid}: no handler registered for kind={kind}. "
                f"Add one to handlers, e.g. handlers['{kind}']=my_handler"
            )
        if sid in samples:
            raise GatherSamplesError(f"Duplicate sample_id '{sid}' after processing. Check your filenames and the aggressive_ids setting.")
        samples[sid] = handler(sid, items)

    return samples


def gather_samples_config(
    inputs: Iterable[Union[str, Path]],
    *,
    max_depth: int = 4,
    aggressive_ids: bool = False,
    exclude_name_regex: Optional[str] = r"^(Undetermined|NTC|PTC)",
    exclude_files: Optional[List[Path]] = None,
) -> Dict:
    """
    Scan inputs for sequence files, build config dict, and return it.
    """
    seqfiles = scan_sequence_files(
        inputs,
        max_depth=max_depth,
        aggressive_ids=aggressive_ids,
        exclude_name_regex=exclude_name_regex,
        exclude_files=exclude_files,
    )
    cfg = build_samples_config(
        seqfiles,
    )
    return cfg

