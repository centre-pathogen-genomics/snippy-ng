import pytest

from tests.integration.simulation import VariantRecord


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
        "strict_region": ("Wildtype", 1, 180),
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
