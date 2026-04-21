import subprocess
import sys
from pathlib import Path
import gzip

import pytest

from snippy_ng.context import Context
from snippy_ng.pipelines.core import CorePipelineBuilder
from tests.integration.simulation import (
    DEFAULT_REFERENCE,
    PROJECT_ROOT,
    SimulationRequest,
    VariantRecord,
    materialize_scenario,
)


pytestmark = pytest.mark.integration_sim


SCENARIOS = [
    {
        "name": "short_snp",
        "input_type": "short",
        "variants": (
            VariantRecord("Wildtype", 20, "A", "C"),
        ),
        "untouched_regions": (("Wildtype", 200, 260),),
        "strict_region": None,
    },
    {
        "name": "short_indel",
        "input_type": "short",
        "variants": (
            VariantRecord("Wildtype", 70, "CA", "C"),
        ),
        "untouched_regions": (("Wildtype", 200, 260),),
        "strict_region": ("Wildtype", 50, 180),
    },
    {
        "name": "long_mixed",
        "input_type": "long",
        "variants": (
            VariantRecord("Wildtype", 120, "A", "C"),
            VariantRecord("Wildtype", 123, "GC", "G"),
        ),
        "untouched_regions": (("Wildtype", 260, 340),),
        "strict_region": None,
    },
    {
        "name": "asm_mixed",
        "input_type": "asm",
        "variants": (
            VariantRecord("Wildtype", 150, "A", "G"),
            VariantRecord("Wildtype", 169, "GT", "G"),
        ),
        "untouched_regions": (("Wildtype", 260, 340),),
        "strict_region": ("Wildtype", 1, 220),
    },
    {
        "name": "negative_region",
        "input_type": "short",
        "variants": (
            VariantRecord("Wildtype", 40, "C", "A"),
        ),
        "untouched_regions": (("Wildtype", 200, 260),),
        "strict_region": None,
    },
    {
        "name": "asm_parametrized",
        "input_type": "asm",
        "variants": (
            VariantRecord("Wildtype", 170, "T", "A"),
            VariantRecord("Wildtype", 190, "T", "G"),
        ),
        "untouched_regions": (("Wildtype", 220, 260),),
        "strict_region": None,
    },
]


@pytest.mark.parametrize(
    ("name", "input_type", "variants", "untouched_regions", "strict_region"),
    [
        (
            scenario["name"],
            scenario["input_type"],
            scenario["variants"],
            scenario["untouched_regions"],
            scenario["strict_region"],
        )
        for scenario in SCENARIOS
    ],
    ids=[scenario["name"] for scenario in SCENARIOS],
)
def test_simulated_pipeline_scenarios(
    simulated_dataset,
    name,
    input_type,
    variants,
    untouched_regions,
    strict_region,
):
    dataset = simulated_dataset(
        variants,
        input_type,
        name=name,
        untouched_regions=untouched_regions,
    )

    for variant in variants:
        dataset.assert_variant_present(variant.chrom, variant.pos, variant.ref, variant.alt)

    for chrom, start, end in dataset.request.untouched_regions:
        dataset.assert_no_variants_in_region(chrom, start, end)

    if strict_region is not None:
        chrom, start, end = strict_region
        dataset.assert_no_unexpected_calls_in_region(
            chrom,
            start,
            end,
            expected=variants,
        )


def test_assembly_pipeline_handles_whole_contig_deletion(tmp_path):
    reference = tmp_path / "reference.fasta"
    assembly = tmp_path / "assembly.fasta"
    outdir = tmp_path / "whole-contig-delete"

    missing_contig = "GATTACA" * 30
    reference.write_text(
        ">contig_1\n"
        + ("ACGT" * 50)
        + "\n>contig_10\n"
        + missing_contig
        + "\n"
    )
    assembly.write_text(">contig_1\n" + ("ACGT" * 50) + "\n")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snippy_ng",
            "asm",
            "--reference",
            str(reference),
            "--assembly",
            str(assembly),
            "--outdir",
            str(outdir),
            "--prefix",
            "snippy",
            "--skip-check",
            "--cpus",
            "1",
            "--ram",
            "4",
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=PROJECT_ROOT,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    deletion_records = []
    with gzip.open(outdir / "snippy.vcf.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] != "contig_10":
                continue
            deletion_records.append(cols)

    assert deletion_records, "Expected a deletion record for the missing contig"
    assert any(
        record[1] == "1"
        and record[4] == "-"
        and f"END={len(missing_contig)}" in record[7]
        and "SVTYPE=DEL" in record[7]
        and "ZERODEPTH" in record[7]
        for record in deletion_records
    ), deletion_records


def _read_fasta_sequence(path: Path, contig: str) -> str:
    current = None
    chunks = []
    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current == contig:
                    return "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
                continue
            if current == contig:
                chunks.append(line)
    if current == contig:
        return "".join(chunks)
    raise AssertionError(f"Contig {contig!r} not found in {path}")


def test_short_consensus_applies_only_pass_variants(tmp_path, integration_cache_root):
    variant = VariantRecord("Wildtype", 20, "A", "C")
    request = SimulationRequest(
        name="short_pass_only_consensus",
        reference=DEFAULT_REFERENCE,
        truth_variants=(variant,),
        untouched_regions=(("Wildtype", 200, 260),),
    )
    materialized = materialize_scenario(
        request=request,
        input_type="short",
        cache_root=integration_cache_root,
    )
    outdir = tmp_path / "short-pass-only-consensus"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snippy_ng",
            "short",
            "--reference", str(request.reference),
            "--R1", str(materialized.reads_r1),
            "--R2", str(materialized.reads_r2),
            "--outdir", str(outdir),
            "--prefix", "snippy",
            "--skip-check",
            "--cpus", "1",
            "--ram", "4",
            "--min-qual", "999999",
            "--depth-mask", "1",
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=PROJECT_ROOT,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    variant_records = []
    with gzip.open(outdir / "snippy.vcf.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] == variant.chrom and cols[1] == str(variant.pos):
                variant_records.append(cols)

    assert variant_records, "Expected called variant record in output VCF"
    assert any(record[6] != "PASS" for record in variant_records), variant_records

    reference_seq = _read_fasta_sequence(request.reference, variant.chrom)
    consensus_seq = _read_fasta_sequence(outdir / "snippy.pseudo.fna", variant.chrom)
    assert consensus_seq[variant.pos - 1] == reference_seq[variant.pos - 1]


def _read_fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current_name: str | None = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(chunks)
                current_name = line[1:].split()[0]
                chunks = []
                continue
            chunks.append(line)
    if current_name is not None:
        records[current_name] = "".join(chunks)
    return records


def _read_phylip_taxa_count(path: Path) -> int:
    with open(path, "r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()
    assert first_line, f"Expected PHYLIP header in {path}"
    return int(first_line.split()[0])


def test_core_pipeline_builds_alignments_and_distance_matrices(simulated_dataset, tmp_path):
    datasets = [
        simulated_dataset(
            (VariantRecord("Wildtype", 20, "A", "C"),),
            "asm",
            name="core_sample_a",
            untouched_regions=(("Wildtype", 200, 260),),
        ),
        simulated_dataset(
            (VariantRecord("Wildtype", 40, "C", "A"),),
            "asm",
            name="core_sample_b",
            untouched_regions=(("Wildtype", 200, 260),),
        ),
        simulated_dataset(
            (VariantRecord("Wildtype", 120, "A", "G"),),
            "asm",
            name="core_sample_c",
            untouched_regions=(("Wildtype", 200, 260),),
        ),
    ]
    outdir = tmp_path / "core-pipeline"

    pipeline = CorePipelineBuilder(
        snippy_dirs=[dataset.outdir for dataset in datasets],
        reference=DEFAULT_REFERENCE,
        prefix="core",
    ).build()
    pipeline.run(
        Context(
            outdir=outdir,
            skip_check=True,
            cpus=1,
            ram=4,
        )
    )

    full_aln = outdir / "core.full.aln"
    soft_core_aln = outdir / "core.095.aln"
    full_phylip = outdir / "core.full.phylip"
    soft_core_phylip = outdir / "core.095.phylip"

    assert full_aln.exists()
    assert soft_core_aln.exists()
    assert full_phylip.exists()
    assert soft_core_phylip.exists()

    full_records = _read_fasta_records(full_aln)
    soft_core_records = _read_fasta_records(soft_core_aln)
    reference_seq = _read_fasta_sequence(DEFAULT_REFERENCE, "Wildtype")

    assert list(full_records) == ["reference", "core_sample_a-asm", "core_sample_b-asm", "core_sample_c-asm"]
    assert list(soft_core_records) == list(full_records)
    assert all(len(seq) == len(reference_seq) for seq in full_records.values())
    assert len(next(iter(soft_core_records.values()))) == 3
    assert {seq.upper() for seq in soft_core_records.values()} == {"ACA", "CCA", "AAA", "ACG"}
    assert _read_phylip_taxa_count(full_phylip) == 4
    assert _read_phylip_taxa_count(soft_core_phylip) == 4
