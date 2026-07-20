from __future__ import annotations

import gzip
import io
import fnmatch
import os
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Iterable, Literal, Mapping, Optional, TypeAlias, TypedDict, Union

from snippy_ng.exceptions import SnippyError


SeqKind: TypeAlias = Literal["ILL", "ONT", "ASM"]
SampleEntry: TypeAlias = dict[str, object]
PathLike: TypeAlias = Union[str, Path]


class SamplesConfig(TypedDict):
    reference: str | None
    samples: dict[str, SampleEntry]


class GatherSamplesError(SnippyError):
    """Raised when sample discovery or sample config construction fails."""


@dataclass(frozen=True)
class SeqFile:
    """A detected sequence-like file and the sample ID inferred from its path."""

    path: Path
    kind: SeqKind
    sample_id: str
    is_alignment: bool = False


SampleHandler: TypeAlias = Callable[[str, list[SeqFile]], SampleEntry]


BIO_COMPRESSION_RE = re.compile(r"\.(?:gz|xz|zstd)$", flags=re.IGNORECASE)
BIO_EXTENSION_RE = re.compile(
    r"\.(?:fq|fastq|fa|fasta|fna|bam|cram|sam|vcf|bcf|gbk|genbank|gff)$",
    flags=re.IGNORECASE,
)
READ_DIRECTION_SUFFIX_RE = re.compile(
    r"(?:_|^)(?:R?[12])(?:(?:[._-][A-Za-z0-9]+)*)$",
    flags=re.IGNORECASE,
)
READ_PAIR_TOKEN_RE = re.compile(
    r"(?:^|_)(R?[12])(?:_\d+)?(?:[._-][A-Za-z0-9]+)*$",
    flags=re.IGNORECASE,
)
ONT_UUID_RE = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$",
    flags=re.IGNORECASE,
)
ONT_KEY_VALUE_RE = re.compile(r"^\w+=\w+")
FASTQ_HEADER_RE = re.compile(r"^@(\S+)")
ALIGNMENT_SUFFIXES = {".bam", ".cram", ".sam"}
LONG_READ_LENGTH = 1_000
ALIGNMENT_PEEK_RECORDS = 1_000
LONG_READ_PLATFORMS = {"ONT", "OXFORDNANOPORE", "NANOPORE", "PACBIO", "PACBI0"}
SHORT_READ_PLATFORMS = {"ILLUMINA", "BGISEQ", "MGISEQ", "DNBSEQ", "SOLID", "LS454"}
LONG_READ_PRESETS = ("lr:", "map-ont", "map-pb", "lr:hq")
SHORT_READ_PRESETS = ("sr:", "short-read", "illumina")


# ---------------------------------------------------------------------------
# Small formatting and path helpers
# ---------------------------------------------------------------------------


def _absolute(path: PathLike) -> Path:
    """Return an absolute path without resolving symlinks."""

    return Path(path).expanduser().absolute()


def _format_paths(paths: Iterable[Path]) -> str:
    return "\n".join(f"- {path}" for path in paths)    


def _raise_duplicate_pair_error(sample_id: str, read_name: str, first: Path, second: Path) -> None:
    raise GatherSamplesError(
        f"{sample_id} (ILL) has multiple {read_name} candidates!\n\n"
        f"Found:\n{_format_paths([first, second])}"
    )


# ---------------------------------------------------------------------------
# Filename and header parsing
# ---------------------------------------------------------------------------


def strip_bio_suffixes(filename: str) -> str:
    """Remove common bioinformatics suffixes and compression suffixes."""

    without_compression = BIO_COMPRESSION_RE.sub("", filename)
    return BIO_EXTENSION_RE.sub("", without_compression)


def strip_read_direction_suffix(sample_id: str) -> str:
    """Remove a terminal Illumina read-direction token from a sample ID."""

    return READ_DIRECTION_SUFFIX_RE.sub("", sample_id)


def read_pair_token(path: Path) -> Optional[Literal["R1", "R2"]]:
    """Return R1/R2 when the filename ends with a recognizable pair token."""

    stem = strip_bio_suffixes(path.name)
    match = READ_PAIR_TOKEN_RE.search(stem)
    if match is None:
        return None

    token = match.group(1).upper()
    return "R1" if token in {"1", "R1"} else "R2"


def _open_text_maybe_gzip(path: Path) -> io.TextIOBase:
    """Open a text stream, handling gzip by extension."""

    if path.suffix.lower() == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def detect_seq_kind(path: PathLike) -> Optional[SeqKind]:
    """
    Detect the sequence kind from the first line of a file.

    Returns
    -------
    "ASM"
        FASTA-like file beginning with ``>``.
    "ONT"
        FASTQ-like file whose read ID looks Oxford Nanopore-like.
    "ILL"
        FASTQ-like file not classified as ONT.
    None
        Missing, empty, unreadable, or unrecognized file.
    """

    file_path = Path(path)
    if not _is_non_empty_file(file_path):
        return None

    if _is_alignment_path(file_path):
        return detect_alignment_seq_kind(file_path)

    try:
        with _open_text_maybe_gzip(file_path) as handle:
            first_line = handle.readline()
    except Exception:
        return None

    if first_line.startswith(">"):
        return "ASM"

    if not first_line.startswith("@"):
        return None

    match = FASTQ_HEADER_RE.match(first_line)
    if match is None:
        return None

    read_id = match.group(1)
    return "ONT" if _looks_like_ont_read_id(read_id) else "ILL"


def detect_alignment_seq_kind(path: PathLike) -> Optional[SeqKind]:
    """Classify an alignment using its platform header when available.

    ``@RG`` ``PL`` tags identify common sequencing platforms without reading
    alignment records. Files without a recognized platform tag fall back to
    peeking at primary read lengths for compatibility with older alignments.
    """

    file_path = Path(path)
    try:
        header_kind = _alignment_header_kind(file_path)
        if header_kind is not None:
            return header_kind

        first_record = _first_alignment_record(file_path)
        if first_record is None or _is_unmapped_alignment(first_record):
            # Skip uBAMs
            return None

        records = _alignment_records(file_path)
        lengths: list[int] = []
        for record in records:
            fields = record.rstrip("\n").split("\t")
            if len(fields) < 11 or fields[9] == "*":
                continue
            try:
                flag = int(fields[1])
                if flag & 0x904:  # unmapped, secondary, supplementary
                    continue
                lengths.append(len(fields[9]))
            except (ValueError, IndexError):
                continue
            if len(lengths) >= ALIGNMENT_PEEK_RECORDS:
                break
    except (OSError, subprocess.SubprocessError):
        return None

    if not lengths:
        return None

    lengths.sort()
    median = lengths[len(lengths) // 2]
    return "ONT" if median >= LONG_READ_LENGTH else "ILL"


def _first_alignment_record(path: Path) -> Optional[str]:
    """Read one alignment record without materializing the whole file."""

    if path.suffix.lower() == ".sam" or path.name.lower().endswith(".sam.gz"):
        with _open_text_maybe_gzip(path) as handle:
            for line in handle:
                if not line.startswith("@"):
                    return line
        return None

    process = subprocess.Popen(
        ["samtools", "view", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        if process.stdout is None:
            return None
        return process.stdout.readline() or None
    finally:
        process.terminate()
        process.wait()


def _is_unmapped_alignment(record: str) -> bool:
    fields = record.rstrip("\n").split("\t")
    if len(fields) < 6:
        return True
    try:
        flag = int(fields[1])
    except ValueError:
        return True
    return bool(flag & 0x4) or fields[5] == "*"


def _alignment_header_kind(path: Path) -> Optional[SeqKind]:
    """Return the sequencing kind encoded by ``@RG`` ``PL`` tags."""

    if path.suffix.lower() == ".sam" or path.name.lower().endswith(".sam.gz"):
        header_lines = list(_sam_header_lines(path))
    else:
        result = subprocess.run(
            ["samtools", "view", "-H", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        header_lines = result.stdout.splitlines()

    header_text = "\n".join(header_lines).lower()
    if any(preset in header_text for preset in LONG_READ_PRESETS):
        return "ONT"
    if any(preset in header_text for preset in SHORT_READ_PRESETS):
        return "ILL"

    platforms = {
        field[3:].upper()
        for line in header_lines
        if line.startswith("@RG")
        for field in line.rstrip("\n").split("\t")
        if field.upper().startswith("PL:")
    }
    if platforms & LONG_READ_PLATFORMS:
        return "ONT"
    if platforms & SHORT_READ_PLATFORMS:
        return "ILL"
    return None


def _sam_header_lines(path: Path) -> Iterable[str]:
    with _open_text_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("@"):
                yield line
            else:
                break


def _alignment_records(path: Path) -> Iterable[str]:
    if path.suffix.lower() == ".sam" or path.name.lower().endswith(".sam.gz"):
        with _open_text_maybe_gzip(path) as handle:
            yield from (line for line in handle if not line.startswith("@"))
        return

    result = subprocess.run(
        ["samtools", "view", "-F", "0x904", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    yield from result.stdout.splitlines()


def _is_alignment_path(path: Path) -> bool:
    name = path.name.lower()
    return path.suffix.lower() in ALIGNMENT_SUFFIXES or name.endswith(".sam.gz")


def _is_non_empty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _looks_like_ont_read_id(read_id: str) -> bool:
    return bool(ONT_UUID_RE.match(read_id) or ONT_KEY_VALUE_RE.match(read_id))


# ---------------------------------------------------------------------------
# Input expansion and reference exclusion
# ---------------------------------------------------------------------------


def iter_candidate_files(inputs: Iterable[PathLike], *, max_depth: int) -> Iterable[Path]:
    """Yield files from input files or directories, respecting max_depth."""

    for input_path in inputs:
        root = _absolute(input_path)

        if root.is_file():
            yield root
            continue

        if root.is_dir():
            yield from _iter_directory_files(root, max_depth=max_depth)


def _iter_directory_files(root: Path, *, max_depth: int) -> Iterable[Path]:
    root_depth = len(root.parts)

    for path in root.rglob("*"):
        candidate = _absolute(path)
        try:
            if not candidate.is_file():
                continue
            if len(candidate.parts) - root_depth > max_depth:
                continue
        except OSError:
            continue

        yield candidate


def collect_reference_exclusions(reference: Optional[PathLike]) -> tuple[Optional[PathLike], list[Path]]:
    """
    Normalize the reference path and return files that should be ignored while scanning.

    A reference can be either a single file or a directory containing metadata.json.
    Directory references exclude every file inside that directory.
    """

    if reference is None:
        return None, []

    from snippy_ng.pipelines.common import is_assembly_accession

    if is_assembly_accession(reference):
        return str(reference), []

    normalized = _absolute(reference)

    if normalized.is_file():
        return normalized, [normalized]

    if normalized.is_dir():
        metadata = normalized / "metadata.json"
        if not metadata.is_file():
            raise GatherSamplesError(
                f"Reference path '{normalized}' is a directory but does not contain metadata.json."
            )
        return normalized, [_absolute(path) for path in normalized.rglob("*") if path.is_file()]

    raise GatherSamplesError(f"Reference path '{normalized}' is neither a file nor a directory.")


# ---------------------------------------------------------------------------
# Sequence file scanning and ID normalization
# ---------------------------------------------------------------------------


def scan_sequence_files(
    inputs: Iterable[PathLike],
    max_depth: int = 4,
    exclude_name_regex: Optional[str] = None,
    exclude_files: Optional[list[Path]] = None,
    reserved_ids: Optional[list[str]] = None,
) -> list[SeqFile]:
    """
    Walk input files/folders and return detected sequence files.

    `exclude_name_regex` is applied to the basename only. Shell-style glob
    patterns such as ``*.fa`` are also accepted.

    If `reserved_ids` is provided, discovered files with the same 
    IDs are disambiguated using the parent folder.
    """

    exclude_re = _compile_exclude_pattern(exclude_name_regex)
    excluded = {_absolute(path) for path in exclude_files or []}

    seqfiles: list[SeqFile] = []
    for path in iter_candidate_files(inputs, max_depth=max_depth):
        if not _should_scan_file(path, exclude_re=exclude_re, excluded=excluded):
            continue

        seqfile = _seqfile_from_path(
            path,
            reserved_ids=reserved_ids,
        )
        if seqfile is not None:
            seqfiles.append(seqfile)

    return normalize_paired_illumina_ids(seqfiles)


def _compile_exclude_pattern(pattern: Optional[str]) -> Optional[re.Pattern[str]]:
    if not pattern:
        return None

    try:
        return re.compile(pattern)
    except re.error:
        if any(character in pattern for character in "*?["):
            return re.compile(fnmatch.translate(pattern))
        raise GatherSamplesError(f"Invalid exclusion pattern {pattern!r}. Use a regular expression or glob such as '*.fa'.")


def _reserved_reference_id(reference: Optional[PathLike]) -> Optional[str]:
    if reference is None:
        return None
    return strip_bio_suffixes(Path(reference).name)


def _should_scan_file(path: Path, *, exclude_re: Optional[re.Pattern[str]], excluded: set[Path]) -> bool:
    if path in excluded:
        return False
    if exclude_re and exclude_re.search(path.name):
        return False
    return True


def _seqfile_from_path(
    path: Path,
    *,
    reserved_ids: Optional[list[str]],
) -> Optional[SeqFile]:
    kind = detect_seq_kind(path)
    if kind is None:
        return None

    sample_id = strip_bio_suffixes(path.name) or path.parent.name
    if reserved_ids and sample_id in reserved_ids:
        sample_id = build_disambiguated_id(sample_id, path.parent.name)

    return SeqFile(path=path, kind=kind, sample_id=sample_id, is_alignment=_is_alignment_path(path))


def normalize_paired_illumina_ids(seqfiles: list[SeqFile]) -> list[SeqFile]:
    """
    Give paired Illumina reads the same sample ID.

    This turns sample_R1.fastq.gz and sample_R2.fastq.gz into two SeqFile objects
    whose sample_id is ``sample``.
    """

    pairs_by_path: dict[Path, str] = {}

    for (_parent, base_id), files in _group_illumina_by_parent_and_base(seqfiles).items():
        r1, r2 = find_r1_r2([seqfile.path for seqfile in files])
        if r1 is None or r2 is None:
            continue
        pairs_by_path[r1] = base_id
        pairs_by_path[r2] = base_id

    return [
        replace(seqfile, sample_id=pairs_by_path[seqfile.path])
        if seqfile.path in pairs_by_path
        else seqfile
        for seqfile in seqfiles
    ]


def _group_illumina_by_parent_and_base(seqfiles: list[SeqFile]) -> dict[tuple[Path, str], list[SeqFile]]:
    grouped: dict[tuple[Path, str], list[SeqFile]] = defaultdict(list)

    for seqfile in seqfiles:
        if seqfile.kind != "ILL":
            continue
        base_id = strip_read_direction_suffix(seqfile.sample_id) or seqfile.path.parent.name
        grouped[(_absolute(seqfile.path.parent), base_id)].append(seqfile)

    return dict(grouped)


# ---------------------------------------------------------------------------
# Illumina pairing
# ---------------------------------------------------------------------------


def find_r1_r2(files: list[Path]) -> tuple[Optional[Path], Optional[Path]]:
    """Find the R1 and R2 paths in a list of candidate files."""

    r1: Optional[Path] = None
    r2: Optional[Path] = None

    for path in files:
        pair = read_pair_token(path)

        if pair == "R1":
            if r1 is not None:
                _raise_duplicate_pair_error(path.name, "R1", r1, path)
            r1 = path
        elif pair == "R2":
            if r2 is not None:
                _raise_duplicate_pair_error(path.name, "R2", r2, path)
            r2 = path

    return r1, r2


# ---------------------------------------------------------------------------
# Duplicate sample ID resolution
# ---------------------------------------------------------------------------


def resolve_duplicate_sample_ids(seqfiles: list[SeqFile]) -> list[SeqFile]:
    """Resolve duplicate sample IDs using parent directory labels or filenames."""

    resolved: list[SeqFile] = []

    for sample_id, files_for_sample in group_by_sample_id(seqfiles).items():
        by_parent = group_by_parent(files_for_sample)

        if len(by_parent) == 1:
            resolved.extend(_resolve_single_parent_group(files_for_sample))
            continue

        parent_labels = build_parent_labels(list(by_parent.keys()))
        for files_in_parent, parent_label in zip(by_parent.values(), parent_labels):
            resolved.extend(_resolve_multi_parent_group(sample_id, files_in_parent, parent_label))

    return resolved


def _resolve_single_parent_group(seqfiles: list[SeqFile]) -> list[SeqFile]:
    if needs_filename_fallback(seqfiles):
        return [replace(seqfile, sample_id=seqfile.path.name) for seqfile in seqfiles]
    return seqfiles


def _resolve_multi_parent_group(
    sample_id: str,
    seqfiles: list[SeqFile],
    parent_label: str,
) -> list[SeqFile]:
    if needs_filename_fallback(seqfiles):
        return [replace(seqfile, sample_id=seqfile.path.name) for seqfile in seqfiles]

    new_id = build_disambiguated_id(sample_id, parent_label)
    return [replace(seqfile, sample_id=new_id) for seqfile in seqfiles]


def group_by_sample_id(seqfiles: Iterable[SeqFile]) -> dict[str, list[SeqFile]]:
    grouped: dict[str, list[SeqFile]] = defaultdict(list)
    for seqfile in seqfiles:
        grouped[seqfile.sample_id].append(seqfile)
    return dict(grouped)


def group_by_parent(seqfiles: Iterable[SeqFile]) -> dict[Path, list[SeqFile]]:
    grouped: dict[Path, list[SeqFile]] = defaultdict(list)
    for seqfile in seqfiles:
        grouped[_absolute(seqfile.path.parent)].append(seqfile)
    return dict(grouped)


def needs_filename_fallback(seqfiles: list[SeqFile]) -> bool:
    """
    True when a group has multiple files but is not a legitimate Illumina group.
    """

    return len(seqfiles) > 1 and not is_expected_multi_file_sample(seqfiles)


def is_expected_multi_file_sample(seqfiles: list[SeqFile]) -> bool:
    """Illumina is the only sample kind that can legitimately have multiple files."""

    if not seqfiles:
        return False
    return {seqfile.kind for seqfile in seqfiles} == {"ILL"}


def build_parent_labels(parents: list[Path]) -> list[str]:
    """Build stable, minimally qualified labels from parent directories."""

    if not parents:
        return []

    absolute_parents = [_absolute(parent) for parent in parents]
    common_parent = Path(os.path.commonpath([str(parent) for parent in absolute_parents]))

    labels: list[str] = []
    for parent in absolute_parents:
        try:
            relative_parts = parent.relative_to(common_parent).parts
        except ValueError:
            relative_parts = parent.parts

        labels.append("-".join(relative_parts) if relative_parts else parent.name)

    return labels


def build_disambiguated_id(sample_id: str, parent_label: str) -> str:
    if not sample_id:
        return parent_label
    if parent_label == sample_id or parent_label.endswith(f"-{sample_id}"):
        return parent_label
    return f"{parent_label}-{sample_id}"


# ---------------------------------------------------------------------------
# Sample entry builders
# ---------------------------------------------------------------------------


def handle_ILL(sample_id: str, seqfiles: list[SeqFile]) -> SampleEntry:
    """Build a config entry for Illumina short-read data."""

    if any(seqfile.is_alignment for seqfile in seqfiles) and len(seqfiles) != 1:
        raise GatherSamplesError(f"{sample_id} (ILL) sample has multiple alignment files!")
    alignment = _single_alignment(seqfiles)
    if alignment is not None:
        return {"type": "short", "bam": str(alignment)}

    paths = sorted(seqfile.path for seqfile in seqfiles)
    r1, r2 = find_r1_r2(paths)

    if r1 is None or r2 is None:
        r1, r2 = _fallback_illumina_pair(sample_id, paths)

    return {
        "type": "short",
        "left": str(r1),
        "right": str(r2) if r2 is not None else None,
    }


def _fallback_illumina_pair(sample_id: str, paths: list[Path]) -> tuple[Path, Optional[Path]]:
    """
    Handle Illumina files whose names do not encode R1/R2.

    One file is treated as single-end. Two files are paired in sort order.
    More than two files is ambiguous.
    """

    if len(paths) == 1:
        return paths[0], None
    if len(paths) == 2:
        return paths[0], paths[1]

    raise GatherSamplesError(
        f"{sample_id}: ILL sample needs 2 files (R1/R2)!\n\n"
        f"Found:\n{_format_paths(paths)}"
    )


def handle_ONT(sample_id: str, seqfiles: list[SeqFile]) -> SampleEntry:
    """Build a config entry for Oxford Nanopore long-read data."""

    if any(seqfile.is_alignment for seqfile in seqfiles) and len(seqfiles) != 1:
        raise GatherSamplesError(f"{sample_id} (ONT) sample has multiple alignment files!")
    alignment = _single_alignment(seqfiles)
    if alignment is not None:
        return {
            "type": "long",
            "bam": str(alignment),
            "caller": "clair3",
            "model": None,
        }

    path = _require_single_file(sample_id, "ONT", seqfiles)
    return {
        "type": "long",
        "reads": str(path),
        "caller": "clair3",
        "model": None,
    }


def _single_alignment(seqfiles: list[SeqFile]) -> Optional[Path]:
    alignments = [seqfile.path for seqfile in seqfiles if seqfile.is_alignment]
    return alignments[0] if len(alignments) == 1 and len(seqfiles) == 1 else None


def handle_ASM(sample_id: str, seqfiles: list[SeqFile]) -> SampleEntry:
    """Build a config entry for assembly data."""

    path = _require_single_file(sample_id, "ASM", seqfiles)
    return {
        "type": "asm",
        "assembly": str(path),
    }


def _require_single_file(sample_id: str, kind: SeqKind, seqfiles: list[SeqFile]) -> Path:
    paths = sorted(seqfile.path for seqfile in seqfiles)

    if not paths:
        raise GatherSamplesError(f"{sample_id} ({kind}) sample has no files")
    if len(paths) > 1:
        raise GatherSamplesError(
        f"{sample_id} ({kind}) sample has multiple files!\n\n"
        f"Found:\n{_format_paths(paths)}"
    )

    return paths[0]


SAMPLE_HANDLERS: Mapping[SeqKind, SampleHandler] = {
    "ILL": handle_ILL,
    "ONT": handle_ONT,
    "ASM": handle_ASM,
}


# ---------------------------------------------------------------------------
# Config construction
# ---------------------------------------------------------------------------


def build_samples_config(seqfiles: list[SeqFile]) -> dict[str, SampleEntry]:
    """Build the ``samples`` section of the config."""

    samples: dict[str, SampleEntry] = {}
    deduped = resolve_duplicate_sample_ids(seqfiles)

    for sample_id, files_for_sample in sorted(group_by_sample_id(deduped).items()):
        kind = require_single_kind(sample_id, files_for_sample)
        handler = SAMPLE_HANDLERS.get(kind)

        if handler is None:
            raise KeyError(
                f"{sample_id}: no handler registered for kind={kind}. "
                f"Add one to SAMPLE_HANDLERS, e.g. SAMPLE_HANDLERS['{kind}'] = my_handler"
            )

        if sample_id in samples:
            raise GatherSamplesError(
                f"Duplicate sample_id '{sample_id}' after processing. "
                "Please check your filenames and ensure they are unique."
            )

        samples[sample_id] = handler(sample_id, files_for_sample)

    return samples


def require_single_kind(sample_id: str, seqfiles: list[SeqFile]) -> SeqKind:
    kinds = {seqfile.kind for seqfile in seqfiles}
    if len(kinds) == 1:
        return next(iter(kinds))

    offending = "\n".join(f"- {seqfile.path} ({seqfile.kind})" for seqfile in seqfiles)
    raise GatherSamplesError(
        f"Multiple kinds detected for sample_id '{sample_id}'.\n\n"
        "Filename parsing must not be ambiguous. "
        "Please rename files or consider adjusting exclude_name_regex settings.\n\n"
        f"Found:\n{offending}"
    )


def apply_sample_defaults(
    samples: dict[str, SampleEntry],
    defaults: Optional[Mapping[str, object]] = None,
) -> dict[str, SampleEntry]:
    """Apply default values to each sample when the key is missing or None."""

    if not defaults:
        return samples

    for sample_data in samples.values():
        for key, value in defaults.items():
            if sample_data.get(key) is None:
                sample_data[key] = value

    return samples


def gather(
    inputs: Iterable[PathLike],
    *,
    max_depth: int = 4,
    exclude_name_regex: Optional[str] = r"^(Undetermined|NTC|PTC)",
    exclude_files: Optional[list[Path]] = None,
    reference: Optional[PathLike] = None,
    defaults: Optional[Mapping[str, object]] = None,
) -> SamplesConfig:
    """Scan inputs, build sample config, apply defaults, and return the final config."""

    normalized_reference, reference_exclusions = collect_reference_exclusions(reference)
    excluded = [_absolute(path) for path in exclude_files or []] + reference_exclusions
    reference_id = _reserved_reference_id(reference)

    reserved_ids = [reference_id] if reference_id is not None else None
    seqfiles = scan_sequence_files(
        inputs,
        max_depth=max_depth,
        exclude_name_regex=exclude_name_regex,
        exclude_files=excluded,
        reserved_ids=reserved_ids,
    )

    samples = build_samples_config(seqfiles)
    apply_sample_defaults(samples, defaults)

    return {
        "reference": str(normalized_reference) if normalized_reference else None,
        "samples": samples,
    }
