import subprocess
import sys

import pytest

from tests.integration.simulation import PROJECT_ROOT, VariantRecord


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
    with open(outdir / "snippy.vcf", "r") as handle:
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
