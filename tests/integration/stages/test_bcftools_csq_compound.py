from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from snippy_ng.context import Context
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller


pytestmark = pytest.mark.integration_sim


MINI_GFF = """##gff-version 3
mini\tsnippy-ng\tgene\t1\t2217\t.\t-\t.\tID=gene:BFV63_RS17260;biotype=protein_coding;Name=hypF
mini\tsnippy-ng\ttranscript\t1\t2217\t.\t-\t.\tID=transcript:BFV63_RS17260_transcript_1;Parent=gene:BFV63_RS17260;biotype=protein_coding
mini\tsnippy-ng\tCDS\t1\t2217\t.\t-\t0\tParent=transcript:BFV63_RS17260_transcript_1
"""

MINI_VCF = """##fileformat=VCFv4.2
##contig=<ID=mini,length=2217>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsnippy
mini\t2214\t.\tT\tG\t60\tPASS\t.\tGT\t1/1
mini\t2216\t.\tA\tAT\t60\tPASS\t.\tGT\t1/1
"""

UNANNOTATED_CONTIG_GFF = """##gff-version 3
annotated\treproducer\tgene\t1\t9\t.\t+\t.\tID=gene1
annotated\treproducer\tmRNA\t1\t9\t.\t+\t.\tID=transcript1;Parent=gene1
annotated\treproducer\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=transcript1
"""

UNANNOTATED_CONTIG_VCF = """##fileformat=VCFv4.2
##contig=<ID=annotated,length=9>
##contig=<ID=contig00032,length=12>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
contig00032\t2\t.\tC\tT\t100\tPASS\t.
"""


def _parse_info_field(info: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for field in info.split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
        else:
            values[field] = ""
    return values


def _annotated_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            records.append(
                {
                    "chrom": cols[0],
                    "pos": cols[1],
                    "ref": cols[3],
                    "alt": cols[4],
                    "filter": cols[6],
                    "info": cols[7],
                    "format": cols[8] if len(cols) > 8 else "",
                    "sample": cols[9] if len(cols) > 9 else "",
                }
            )
    return records


def _build_reference_sequence(length: int = 2217, seed: int = 7) -> str:
    rng = random.Random(seed)
    prefix = "".join(rng.choice("ACGT") for _ in range(length - 6))
    return prefix + "GCTCAT"


def _write_compound_bug_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    reference = tmp_path / "mini.fa"
    features = tmp_path / "mini.gff"
    variants = tmp_path / "mini.vcf"

    reference.write_text(f">mini\n{_build_reference_sequence()}\n", encoding="utf-8")
    features.write_text(MINI_GFF, encoding="utf-8")
    variants.write_text(MINI_VCF, encoding="utf-8")

    return reference, features, variants


def _write_unannotated_contig_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    reference = tmp_path / "unannotated-contig.fa"
    features = tmp_path / "unannotated-contig.gff"
    variants = tmp_path / "unannotated-contig.vcf"

    reference.write_text(">annotated\nATGAAATAA\n>contig00032\nACGTACGTACGT\n", encoding="utf-8")
    features.write_text(UNANNOTATED_CONTIG_GFF, encoding="utf-8")
    variants.write_text(UNANNOTATED_CONTIG_VCF, encoding="utf-8")

    return reference, features, variants


def test_bcftools_csq_compound_reverse_strand_bug_reproducer_is_annotated(tmp_path: Path):
    """
    There is a known bug in bcftools csq that can cause incorrect annotations for compound variants on the reverse strand when not using --local-csq mode
    https://github.com/samtools/bcftools/issues/2543
    """
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools required for this consequence regression test")

    reference, features, variants = _write_compound_bug_fixture(tmp_path)
    stage = BcftoolsConsequencesCaller(
        prefix=str(tmp_path / "compound-csq"),
        reference=reference,
        variants=variants,
        features=features,
    )

    stage.run(Context(cpus=1, quiet=True))
    stage.run_tests()

    records = _annotated_records(stage.output.annotated_vcf)

    assert len(records) == 2
    assert [(record["chrom"], record["pos"], record["ref"], record["alt"]) for record in records] == [
        ("mini", "2214", "T", "G"),
        ("mini", "2216", "A", "AT"),
    ]
    assert all(record["filter"] == "PASS" for record in records)
    assert all(record["format"] == "GT:BCSQ" for record in records)
    assert all(record["sample"] == "1/1:3" for record in records)

    annotations = [_parse_info_field(record["info"])["BCSQ"] for record in records]
    assert annotations == [
        "missense|hypF|BFV63_RS17260_transcript_1|protein_coding|-|2S>2R|2214T>G",
        "start_lost|hypF|BFV63_RS17260_transcript_1|protein_coding|-",
    ]


def test_bcftools_csq_unannotated_contig_variant_passes_through(tmp_path: Path):
    """
    bcftools 1.23 exits unless --force is used when a VCF has variants on a
    valid reference contig with no GFF feature records.
    https://github.com/samtools/bcftools/issues/2584
    """
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools required for this consequence regression test")

    reference, features, variants = _write_unannotated_contig_fixture(tmp_path)
    stage = BcftoolsConsequencesCaller(
        prefix=str(tmp_path / "unannotated-contig-csq"),
        reference=reference,
        variants=variants,
        features=features,
    )

    stage.run(Context(cpus=1, quiet=True))
    stage.run_tests()

    records = _annotated_records(stage.output.annotated_vcf)

    assert len(records) == 1
    assert records[0]["chrom"] == "contig00032"
    assert records[0]["pos"] == "2"
    assert records[0]["ref"] == "C"
    assert records[0]["alt"] == "T"
    assert records[0]["filter"] == "PASS"
    assert records[0]["info"] == "."
