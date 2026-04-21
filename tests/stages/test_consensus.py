import gzip
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import pytest

from snippy_ng.context import Context
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
DISTANT_REFERENCE = DATA_DIR / "distant.fa"
SNIPPY_PASS_VCF = DATA_DIR / "snippy.pass.vcf"

IUPAC_CODES = {
    frozenset(("A", "C")): "M",
    frozenset(("A", "G")): "R",
    frozenset(("A", "T")): "W",
    frozenset(("C", "G")): "S",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "T")): "K",
}


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    filter: str
    info: str
    gt: str | None


def _load_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_name: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
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


def _parse_vcf_records(path: Path) -> list[VcfRecord]:
    records: list[VcfRecord] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            fmt = cols[8].split(":") if len(cols) > 8 else []
            sample = cols[9].split(":") if len(cols) > 9 else []
            gt = sample[fmt.index("GT")] if "GT" in fmt and len(sample) > fmt.index("GT") else None
            records.append(
                VcfRecord(
                    chrom=cols[0],
                    pos=int(cols[1]),
                    ref=cols[3].upper(),
                    alt=cols[4].split(",")[0].upper(),
                    filter=cols[6],
                    info=cols[7],
                    gt=gt,
                )
            )
    return records


def _is_homozygous_alt(gt: str | None) -> bool:
    if gt is None:
        return False
    alleles = re.split(r"[/|]", gt)
    return bool(alleles) and all(allele == "1" for allele in alleles)


def _is_heterozygous_ref_alt(gt: str | None) -> bool:
    if gt is None:
        return False
    alleles = re.split(r"[/|]", gt)
    return sorted(alleles) == ["0", "1"]


def _info_value(info: str, key: str) -> int | None:
    for field in info.split(";"):
        if field.startswith(f"{key}="):
            return int(field.split("=", 1)[1])
    return None


def _iupac_consensus(ref: str, alt: str) -> str | None:
    chars: list[str] = []
    for ref_base, alt_base in zip(ref, alt):
        if ref_base == alt_base:
            chars.append(ref_base)
        elif ref_base in "ACGT" and alt_base in "ACGT":
            chars.append(IUPAC_CODES[frozenset((ref_base, alt_base))])
        else:
            return None
    return "".join(chars)


def _expected_consensus_span(record: VcfRecord) -> tuple[int, str] | None:
    """Return zero-based start and expected fixed-length consensus sequence."""
    if record.filter != "PASS":
        return None
    if len(record.alt) > len(record.ref) and record.alt != "<DEL>":
        return None

    if record.alt == "<DEL>":
        if not _is_homozygous_alt(record.gt):
            return None
        end = _info_value(record.info, "END")
        if end is None:
            return None
        # Symbolic DEL records use a left anchor: POS is retained, POS+1..END is deleted.
        return record.pos, "-" * (end - record.pos)

    if _is_homozygous_alt(record.gt):
        return record.pos - 1, record.alt + ("-" * (len(record.ref) - len(record.alt)))

    if _is_heterozygous_ref_alt(record.gt) and len(record.ref) == len(record.alt):
        expected = _iupac_consensus(record.ref, record.alt)
        if expected is not None:
            return record.pos - 1, expected

    return None


def test_bcftools_consensus_applies_all_snippy_variants(tmp_path: Path):
    missing_commands = [command for command in ("bcftools", "bgzip") if shutil.which(command) is None]
    if missing_commands:
        pytest.skip(f"{', '.join(missing_commands)} required for this consensus regression test")
    if not SNIPPY_PASS_VCF.exists():
        pytest.skip(f"{SNIPPY_PASS_VCF} is required for this consensus regression test")

    reference_sequences = _load_fasta(DISTANT_REFERENCE)
    pass_vcf = tmp_path / SNIPPY_PASS_VCF.name
    shutil.copyfile(SNIPPY_PASS_VCF, pass_vcf)

    stage = BcftoolsPseudoAlignment(
        prefix=str(tmp_path / "snippy"),
        reference=DISTANT_REFERENCE,
        vcf=pass_vcf,
        ref_metadata=ReferenceMetadata(
            total_length=sum(len(seq) for seq in reference_sequences.values()),
            num_sequences=len(reference_sequences),
        ),
    )
    stage.run(Context(quiet=True))
    stage.run_tests()

    consensus_sequences = _load_fasta(stage.output.fasta)
    mismatches: list[str] = []
    checked = 0
    for record in _parse_vcf_records(pass_vcf):
        expected_span = _expected_consensus_span(record)
        if expected_span is None:
            continue
        start, expected = expected_span
        observed = consensus_sequences[record.chrom][start : start + len(expected)]
        checked += 1
        if observed != expected:
            mismatches.append(
                f"{record.chrom}:{record.pos} {record.ref}>{record.alt} GT={record.gt} "
                f"expected {expected!r}, observed {observed!r}"
            )

    assert checked > 0
    assert not mismatches, (
        f"{len(mismatches)} of {checked} representable PASS non-insertion variants "
        "were not applied correctly:\n" + "\n".join(mismatches[:20])
    )
