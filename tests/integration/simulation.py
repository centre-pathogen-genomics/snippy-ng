from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CACHE_ROOT = PROJECT_ROOT / ".cache" / "integration-sim"
SCENARIO_DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_REFERENCE = SCENARIO_DATA_DIR / "reference.fasta"


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass(frozen=True)
class SimulationRequest:
    name: str
    reference: Path
    truth_variants: tuple[VariantRecord, ...]
    untouched_regions: tuple[tuple[str, int, int], ...] = ()
    short_coverage: int = 40
    short_read_length: int = 150
    short_fragment_mean: int = 350
    short_fragment_std: int = 25
    short_seed: int = 7
    long_coverage: int = 20
    long_seed: int = 11
    long_length_mean: int = 900
    long_length_stdev: int = 100
    long_identity_mean: float = 97.0
    long_identity_stdev: float = 1.5


@dataclass(frozen=True)
class MaterializedScenario:
    request: SimulationRequest
    input_type: str
    cache_dir: Path
    mutated_reference: Path
    truth_vcf: Path
    assembly: Path | None = None
    reads_r1: Path | None = None
    reads_r2: Path | None = None
    long_reads: Path | None = None
    manifest: Path | None = None


@dataclass
class IntegrationDataset:
    request: SimulationRequest
    input_type: str
    reference: Path
    truth_vcf: Path
    mutated_reference: Path
    outdir: Path
    called_vcf: Path
    assembly: Path | None = None
    reads_r1: Path | None = None
    reads_r2: Path | None = None
    long_reads: Path | None = None

    def called_records(self) -> list[VariantRecord]:
        return parse_vcf_records(self.called_vcf)

    def assert_variant_present(self, chrom: str, pos: int, ref: str, alt: str) -> None:
        expected = normalize_variant(VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt))
        called = {normalize_variant(record) for record in self.called_records()}
        assert expected in called, f"Expected variant {expected} not found in {self.called_vcf}"

    def assert_no_variants_in_region(self, chrom: str, start: int, end: int) -> None:
        unexpected = [
            record for record in self.called_records()
            if record.chrom == chrom and start <= record.pos <= end
        ]
        assert not unexpected, f"Unexpected variants in {chrom}:{start}-{end}: {unexpected}"

    def assert_no_unexpected_calls_in_region(
        self,
        chrom: str,
        start: int,
        end: int,
        expected: Iterable[VariantRecord],
    ) -> None:
        expected_set = {
            normalize_variant(record)
            for record in expected
            if record.chrom == chrom and start <= record.pos <= end
        }
        called_set = {
            normalize_variant(record)
            for record in self.called_records()
            if record.chrom == chrom and start <= record.pos <= end
        }
        extra = sorted(called_set - expected_set, key=lambda record: (record.chrom, record.pos, record.ref, record.alt))
        assert not extra, f"Unexpected variants in {chrom}:{start}-{end}: {extra}"


REQUIRED_COMMANDS: dict[str, tuple[str, ...]] = {
    "short": ("art_illumina", "samtools", "bcftools", "freebayes-parallel", "fasta_generate_regions.py", "minimap2", "bedtools"),
    "long": ("badread", "samtools", "bcftools", "freebayes-parallel", "fasta_generate_regions.py", "minimap2", "bedtools", "seqkit"),
    "asm": ("samtools", "bcftools", "minimap2", "paftools.js", "bedtools"),
}


def ensure_commands_available(input_type: str) -> None:
    missing = [command for command in REQUIRED_COMMANDS[input_type] if shutil.which(command) is None]
    if missing:
        commands = ", ".join(missing)
        raise RuntimeError(f"Missing required commands for {input_type} integration test: {commands}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tool_version(command: str) -> str:
    probes = (
        [command, "--version"],
        [command, "version"],
        [command],
    )
    for probe in probes:
        try:
            result = subprocess.run(probe, check=False, capture_output=True, text=True)
        except FileNotFoundError:
            return "missing"
        text = "\n".join(part.strip() for part in (result.stdout, result.stderr) if part.strip())
        if text:
            return text.splitlines()[0].strip()
    return "unknown"


def build_cache_key(request: SimulationRequest, input_type: str) -> str:
    truth_records = [
        {
            "chrom": record.chrom,
            "pos": record.pos,
            "ref": record.ref,
            "alt": record.alt,
        }
        for record in request.truth_variants
    ]
    payload = {
        "name": request.name,
        "input_type": input_type,
        "reference_sha256": sha256_file(request.reference),
        "truth_records": truth_records,
        "untouched_regions": request.untouched_regions,
        "short_coverage": request.short_coverage,
        "short_read_length": request.short_read_length,
        "short_fragment_mean": request.short_fragment_mean,
        "short_fragment_std": request.short_fragment_std,
        "short_seed": request.short_seed,
        "long_coverage": request.long_coverage,
        "long_seed": request.long_seed,
        "long_length_mean": request.long_length_mean,
        "long_length_stdev": request.long_length_stdev,
        "long_identity_mean": request.long_identity_mean,
        "long_identity_stdev": request.long_identity_stdev,
        "tool_versions": {
            command: tool_version(command)
            for command in REQUIRED_COMMANDS.get(input_type, ())
            if command in {"art_illumina", "badread"}
        },
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def parse_vcf_records(vcf_path: Path) -> list[VariantRecord]:
    records: list[VariantRecord] = []
    with open(vcf_path, "r") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            alt = cols[4].split(",")[0]
            if alt in {".", "*"}:
                continue
            records.append(
                VariantRecord(
                    chrom=cols[0],
                    pos=int(cols[1]),
                    ref=cols[3],
                    alt=alt,
                )
            )
    return records


def normalize_variant(record: VariantRecord) -> VariantRecord:
    pos = record.pos
    ref = record.ref
    alt = record.alt
    if alt.startswith("<") and alt.endswith(">"):
        return record

    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return VariantRecord(chrom=record.chrom, pos=pos, ref=ref, alt=alt)


def _load_fasta(reference_fasta: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_name: str | None = None
    chunks: list[str] = []
    with open(reference_fasta, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(chunks).upper()
                current_name = line[1:].split()[0]
                chunks = []
                continue
            chunks.append(line)
    if current_name is not None:
        sequences[current_name] = "".join(chunks).upper()
    return sequences


def _write_fasta(sequences: dict[str, str], output_fasta: Path) -> None:
    with open(output_fasta, "w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx:idx + 80] + "\n")


def write_truth_vcf(records: Iterable[VariantRecord], output_vcf: Path) -> None:
    with open(output_vcf, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for record in records:
            handle.write(
                f"{record.chrom}\t{record.pos}\t.\t{record.ref}\t{record.alt}\t.\tPASS\t.\n"
            )


def apply_truth_vcf(reference_fasta: Path, truth_vcf: Path, output_fasta: Path) -> None:
    sequences = _load_fasta(reference_fasta)
    truth_records = parse_vcf_records(truth_vcf)
    edits_by_chrom: dict[str, list[VariantRecord]] = {}
    for record in truth_records:
        if record.alt.startswith("<") and record.alt.endswith(">"):
            raise ValueError(f"Symbolic alleles are not supported in simulation truth VCFs: {record}")
        edits_by_chrom.setdefault(record.chrom, []).append(record)

    mutated_sequences: dict[str, str] = {}
    for chrom, sequence in sequences.items():
        records = sorted(edits_by_chrom.get(chrom, []), key=lambda record: record.pos)
        cursor = 0
        chunks: list[str] = []
        for record in records:
            start = record.pos - 1
            end = start + len(record.ref)
            if start < cursor:
                raise ValueError(f"Overlapping variants are not supported: {record}")
            if sequence[start:end].upper() != record.ref.upper():
                raise ValueError(
                    f"Reference mismatch for {chrom}:{record.pos}: expected {record.ref}, found {sequence[start:end]}"
                )
            chunks.append(sequence[cursor:start])
            chunks.append(record.alt)
            cursor = end
        chunks.append(sequence[cursor:])
        mutated_sequences[chrom] = "".join(chunks)

    missing_contigs = set(edits_by_chrom) - set(sequences)
    if missing_contigs:
        missing = ", ".join(sorted(missing_contigs))
        raise ValueError(f"Truth VCF references contigs missing from FASTA: {missing}")

    _write_fasta(mutated_sequences, output_fasta)


def run_command(args: list[str], stdout_path: Path | None = None) -> None:
    stdout_handle = None
    try:
        if stdout_path is not None:
            stdout_handle = open(stdout_path, "w")
        subprocess.run(
            args,
            check=True,
            stdout=stdout_handle,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Command failed: {' '.join(args)}\nSTDERR:\n{exc.stderr}"
        ) from exc
    finally:
        if stdout_handle is not None:
            stdout_handle.close()


def simulate_short_reads(request: SimulationRequest, mutated_reference: Path, output_prefix: Path) -> tuple[Path, Path]:
    run_command(
        [
            "art_illumina",
            "-ss", "HS25",
            "-i", str(mutated_reference),
            "-p",
            "-l", str(request.short_read_length),
            "-f", str(request.short_coverage),
            "-m", str(request.short_fragment_mean),
            "-s", str(request.short_fragment_std),
            "-rs", str(request.short_seed),
            "-na",
            "-o", str(output_prefix),
        ]
    )
    return output_prefix.parent / f"{output_prefix.name}1.fq", output_prefix.parent / f"{output_prefix.name}2.fq"


def simulate_long_reads(request: SimulationRequest, mutated_reference: Path, output_fastq: Path) -> Path:
    run_command(
        [
            "badread",
            "simulate",
            "--reference", str(mutated_reference),
            "--quantity", f"{request.long_coverage}x",
            "--length", f"{request.long_length_mean},{request.long_length_stdev}",
            "--identity", f"{request.long_identity_mean},{request.long_identity_mean + request.long_identity_stdev},{request.long_identity_stdev}",
            "--error_model", "random",
            "--qscore_model", "ideal",
            "--seed", str(request.long_seed),
        ],
        stdout_path=output_fastq,
    )
    return output_fastq


def materialize_scenario(
    request: SimulationRequest,
    input_type: str,
    cache_root: Path = DEFAULT_CACHE_ROOT,
) -> MaterializedScenario:
    if input_type not in {"short", "long", "asm"}:
        raise ValueError(f"Unsupported input_type: {input_type}")

    ensure_commands_available(input_type)
    cache_key = build_cache_key(request, input_type)
    cache_dir = cache_root / request.name / input_type / cache_key
    manifest = cache_dir / "manifest.json"
    truth_vcf = cache_dir / f"{request.name}.truth.vcf"
    mutated_reference = cache_dir / f"{request.name}.mutated.fasta"
    assembly = cache_dir / f"{request.name}.assembly.fasta"
    reads_r1 = cache_dir / f"{request.name}.R1.fq"
    reads_r2 = cache_dir / f"{request.name}.R2.fq"
    long_reads = cache_dir / f"{request.name}.long.fastq"

    expected_paths = {
        "short": (truth_vcf, mutated_reference, reads_r1, reads_r2),
        "long": (truth_vcf, mutated_reference, long_reads),
        "asm": (truth_vcf, mutated_reference, assembly),
    }[input_type]

    if manifest.exists() and all(path.exists() for path in expected_paths):
        return MaterializedScenario(
            request=request,
            input_type=input_type,
            cache_dir=cache_dir,
            mutated_reference=mutated_reference,
            truth_vcf=truth_vcf,
            assembly=assembly if input_type == "asm" else None,
            reads_r1=reads_r1 if input_type == "short" else None,
            reads_r2=reads_r2 if input_type == "short" else None,
            long_reads=long_reads if input_type == "long" else None,
            manifest=manifest,
        )

    cache_dir.mkdir(parents=True, exist_ok=True)
    write_truth_vcf(request.truth_variants, truth_vcf)
    apply_truth_vcf(request.reference, truth_vcf, mutated_reference)

    if input_type == "short":
        simulated_r1, simulated_r2 = simulate_short_reads(request, mutated_reference, cache_dir / f"{request.name}.R")
        if simulated_r1 != reads_r1:
            simulated_r1.replace(reads_r1)
        if simulated_r2 != reads_r2:
            simulated_r2.replace(reads_r2)
    elif input_type == "long":
        simulate_long_reads(request, mutated_reference, long_reads)
    else:
        shutil.copyfile(mutated_reference, assembly)

    manifest.write_text(
        json.dumps(
            {
                "name": request.name,
                "input_type": input_type,
                "cache_key": cache_key,
                "reference": str(request.reference),
                "truth_vcf": str(truth_vcf),
            },
            indent=2,
            sort_keys=True,
        )
    )

    return MaterializedScenario(
        request=request,
        input_type=input_type,
        cache_dir=cache_dir,
        mutated_reference=mutated_reference,
        truth_vcf=truth_vcf,
        assembly=assembly if input_type == "asm" else None,
        reads_r1=reads_r1 if input_type == "short" else None,
        reads_r2=reads_r2 if input_type == "short" else None,
        long_reads=long_reads if input_type == "long" else None,
        manifest=manifest,
    )
