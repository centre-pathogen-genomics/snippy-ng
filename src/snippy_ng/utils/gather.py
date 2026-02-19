from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Literal, Optional, Tuple, Union
from snippy_ng.exceptions import SnippyError
import gzip
import io
import os
import re

class GatherSamplesError(SnippyError):
    pass


def _format_file_list(paths: Iterable[Path]) -> str:
    """Format a list of paths as bullet points."""
    return '\n'.join(f"- {str(p)}" for p in paths)


def _raise_multiple_files_error(sample_id: str, sample_type: str, paths: List[Path]) -> None:
    """Raise error for samples that should have exactly one file but have multiple."""
    raise GatherSamplesError(
        f"{sample_id} ({sample_type}) sample has multiple files!\n\n"
        f"Found:\n{_format_file_list(paths)}"
    )


def _raise_duplicate_candidate_error(sample_id: str, file_type: str, first: Path, second: Path) -> None:
    """Raise error when multiple candidates are found for the same file type."""
    raise GatherSamplesError(
        f"{sample_id} (ILL) has multiple {file_type} candidates!\n\n"
        f"Found:\n{_format_file_list([first, second])}"
    )

SeqKind = str  # "ILL" | "ONT" | "ASM" | ...

@dataclass(frozen=True)
class SeqFile:
    path: Path
    kind: SeqKind           # e.g. "ILL", "ONT", "ASM"
    sample_id: str          # inferred ID from filename


def _to_absolute_path(path: Path) -> Path:
    """Return absolute path without resolving symlinks."""
    return Path(path).absolute()


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
    # examples: _R1, _R1_001, _1, _2, R1, R2 (with or without prefix)
    # Match either: underscore prefix OR start of string (for files like R1.fastq)
    name = re.sub(r"(?:_|^)(?:R?[12])(?:_\d+)?$", "", name)

    return name


def scan_sequence_files(
    inputs: Iterable[Union[str, Path]],
    max_depth: int = 4,
    aggressive_ids: bool = False,
    exclude_name_regex: Optional[str] = None,
    exclude_files: Optional[List[Path]] = None,
    reference: Optional[Union[str, Path]] = None,
) -> List[SeqFile]:
    """
    Walk the given files/folders (recursively up to max_depth) and detect sequence files.
    exclude_name_regex is applied to the basename only (not full path).
    If reference is provided, its inferred sample_id is treated as reserved and
    discovered files with the same sample_id are disambiguated by parent folder.
    """
    exclude_re = re.compile(exclude_name_regex) if exclude_name_regex else None
    reserved_reference_id: Optional[str] = None
    if reference is not None:
        reserved_reference_id = guess_sample_id(Path(reference).name, aggressive=aggressive_ids)

    found: List[SeqFile] = []
    for inp in inputs:
        root = Path(inp).expanduser()
        if root.is_file():
            candidates = [root]
        elif root.is_dir():
            # manual depth control
            base_parts = len(_to_absolute_path(root).parts)
            candidates = []
            for p in root.rglob("*"):
                p = _to_absolute_path(p)
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
            if exclude_re and exclude_re.search(str(p)):
                continue
            kind = detect_seq_kind(p)
            if not kind:
                continue
            sid = guess_sample_id(p.name, aggressive=aggressive_ids)
            if not sid:
                sid = p.parent.name
            if reserved_reference_id and sid == reserved_reference_id:
                sid = _build_disambiguated_id(sid, p.parent.name)
            found.append(SeqFile(path=p, kind=kind, sample_id=sid))

    return found


def _find_r1_r2(files: List[Path]) -> Tuple[Optional[Path], Optional[Path]]:
    r1 = None
    r2 = None
    for p in files:
        n = p.name
        if re.search(r"(?:R?1(?:_\d+)?)(?:\.\w+)*$", n):
            if r1 is not None:
                _raise_duplicate_candidate_error(n, "R1", r1, p)
            r1 = p
        elif re.search(r"(?:R?2(?:_\d+)?)(?:\.\w+)*$", n):
            if r2 is not None:
                _raise_duplicate_candidate_error(n, "R2", r2, p)
            r2 = p
    return r1, r2


def handle_ILL(sample_id: str, items: List[SeqFile]) -> Dict:
    # "short": left/right
    paths = sorted([x.path for x in items])
    r1, r2 = _find_r1_r2(paths)

    # Fallback if names don't encode R1/R2: take first two in sort order
    if r1 is None or r2 is None:
        if len(paths) == 1:
            # single file, assume unpaired
            r1 = paths[0]
            r2 = None
        elif len(paths) == 2:
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
        raise GatherSamplesError(f"{sample_id} (ONT) sample has no files")
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
        raise GatherSamplesError(f"{sample_id} (ASM) sample has no files")
    if len(paths) > 1:
        _raise_multiple_files_error(sample_id, "ASM", paths)
    return {
        "type": "asm",
        "assembly": str(paths[0]),
    }


SampleHandler = Callable[[str, List[SeqFile], dict], Dict]


def _is_expected_multi_file_sample(items: List[SeqFile]) -> bool:
    """True when multiple files are expected for a single sample ID."""
    if not items:
        return False
    kinds = {sf.kind for sf in items}
    if len(kinds) != 1:
        return False
    kind = next(iter(kinds))
    # Illumina can be paired-end and legitimately have multiple files.
    return kind == "ILL"


def _build_parent_labels(parents: List[Path]) -> List[str]:
    """
    Build stable, minimally qualified labels for parent directories.
    """
    if not parents:
        return []

    absolute_parents = [_to_absolute_path(p) for p in parents]
    common_parent = Path(os.path.commonpath([str(p) for p in absolute_parents]))
    labels: List[str] = []
    for p in absolute_parents:
        try:
            rel_parts = p.relative_to(common_parent).parts
        except ValueError:
            rel_parts = p.parts
        if rel_parts:
            labels.append("-".join(rel_parts))
        else:
            labels.append(p.name)
    return labels


def _group_by_sample_id(seqfiles: List[SeqFile]) -> Dict[str, List[SeqFile]]:
    grouped: Dict[str, List[SeqFile]] = {}
    for sf in seqfiles:
        grouped.setdefault(sf.sample_id, []).append(sf)
    return grouped


def _group_by_parent(seqfiles: List[SeqFile]) -> Dict[Path, List[SeqFile]]:
    grouped: Dict[Path, List[SeqFile]] = {}
    for sf in seqfiles:
        grouped.setdefault(_to_absolute_path(sf.path.parent), []).append(sf)
    return grouped


def _is_filename_fallback_needed(items: List[SeqFile]) -> bool:
    """
    True for multi-file groups that are not expected ILL pair/single-end groups.
    """
    return len(items) > 1 and not _is_expected_multi_file_sample(items)


def _build_disambiguated_id(sample_id: str, parent_label: str) -> str:
    if not sample_id:
        return parent_label
    if parent_label == sample_id or parent_label.endswith(f"-{sample_id}"):
        return parent_label
    return f"{parent_label}-{sample_id}"


def _resolve_duplicate_ids(seqfiles: List[SeqFile]) -> List[SeqFile]:
    """
    Resolve duplicate sample IDs by incorporating parent directory names.
    """
    # Group by initial sample_id
    grouped = _group_by_sample_id(seqfiles)
    
    deduped_seqfiles: List[SeqFile] = []
    
    for sid, items in grouped.items():
        by_parent = _group_by_parent(items)
        
        if len(by_parent) == 1:
            # Same directory. Keep legitimate multi-file samples as-is.
            if not _is_filename_fallback_needed(items):
                deduped_seqfiles.extend(items)
            else:
                # Disambiguate by keeping filename (with extension), e.g.
                # Sample.fa and Sample.fa.gz remain distinct sample IDs.
                for sf in items:
                    deduped_seqfiles.append(SeqFile(
                        path=sf.path,
                        kind=sf.kind,
                        sample_id=sf.path.name,
                    ))
        else:
            # We have duplicates - need to disambiguate
            parent_dirs = list(by_parent.keys())
            parent_labels = _build_parent_labels(parent_dirs)
            
            # Create new SeqFiles with updated IDs
            for (_, files), parent_label in zip(by_parent.items(), parent_labels):
                # If a parent directory itself has a non-ILL clash (e.g. mixed kinds
                # like JKD6159.fastq.gz + JKD6159.fasta), keep file extensions.
                if _is_filename_fallback_needed(files):
                    for sf in files:
                        deduped_seqfiles.append(SeqFile(
                            path=sf.path,
                            kind=sf.kind,
                            sample_id=sf.path.name,
                        ))
                    continue

                new_id = _build_disambiguated_id(sid, parent_label)
                for sf in files:
                    deduped_seqfiles.append(SeqFile(
                        path=sf.path,
                        kind=sf.kind,
                        sample_id=new_id
                    ))
    
    return deduped_seqfiles


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
    
    # Resolve duplicate IDs first
    seqfiles = _resolve_duplicate_ids(seqfiles)
    
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
    reference: Optional[Union[str, Path]] = None,
) -> Dict[Literal["samples", "reference"], Union[Dict[str, Dict], Optional[str]]]:
    """
    Scan inputs for sequence files, build config dict, and return it.
    """
    normalized_reference = _to_absolute_path(Path(reference)) if reference else None
    normalized_exclude_files: List[Path] = [_to_absolute_path(p) for p in exclude_files] if exclude_files else []
    if normalized_reference:
        normalized_exclude_files.append(normalized_reference)

    seqfiles = scan_sequence_files(
        inputs,
        max_depth=max_depth,
        aggressive_ids=aggressive_ids,
        exclude_name_regex=exclude_name_regex,
        exclude_files=normalized_exclude_files,
        reference=normalized_reference,
    )
    samples = build_samples_config(
        seqfiles,
    )
    return {
        "reference": str(normalized_reference) if normalized_reference else None,
        "samples": samples,
    }

