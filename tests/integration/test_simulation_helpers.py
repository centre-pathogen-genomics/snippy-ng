from pathlib import Path

from tests.integration.simulation import (
    SimulationRequest,
    VariantRecord,
    apply_truth_vcf,
    build_cache_key,
    materialize_scenario,
    normalize_variant,
    parse_vcf_records,
)


def test_apply_truth_vcf_applies_snp_and_deletion(tmp_path: Path):
    reference = tmp_path / "ref.fasta"
    truth_vcf = tmp_path / "truth.vcf"
    output = tmp_path / "mutated.fasta"

    reference.write_text(">ref\nAACCGGTT\n")
    truth_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "ref\t3\t.\tC\tT\t.\tPASS\t.\n"
        "ref\t6\t.\tGT\tG\t.\tPASS\t.\n"
    )

    apply_truth_vcf(reference, truth_vcf, output)

    assert output.read_text() == ">ref\nAATCGGT\n"


def test_normalize_variant_trims_shared_padding():
    assert normalize_variant(VariantRecord("ref", 10, "AC", "AT")) == VariantRecord("ref", 11, "C", "T")


def test_build_cache_key_changes_with_truth_vcf(tmp_path: Path, monkeypatch):
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref\nAACCGGTT\n")
    monkeypatch.setattr("tests.integration.simulation.REQUIRED_COMMANDS", {"asm": ()})

    request_a = SimulationRequest(
        name="demo",
        reference=reference,
        truth_variants=(VariantRecord("ref", 3, "C", "T"),),
    )
    request_b = SimulationRequest(
        name="demo",
        reference=reference,
        truth_variants=(VariantRecord("ref", 4, "C", "A"),),
    )

    assert build_cache_key(request_a, "asm") != build_cache_key(request_b, "asm")


def test_materialize_scenario_skips_regeneration_when_cached(tmp_path: Path, monkeypatch):
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref\nAACCGGTT\n")

    request = SimulationRequest(
        name="demo",
        reference=reference,
        truth_variants=(VariantRecord("ref", 3, "C", "T"),),
    )

    monkeypatch.setattr("tests.integration.simulation.REQUIRED_COMMANDS", {"asm": ()})
    first = materialize_scenario(request, "asm", cache_root=tmp_path / "cache")

    called = {"count": 0}

    def fail_if_called(*args, **kwargs):
        called["count"] += 1
        raise AssertionError("simulation should not run when cache exists")

    monkeypatch.setattr("tests.integration.simulation.apply_truth_vcf", fail_if_called)
    second = materialize_scenario(request, "asm", cache_root=tmp_path / "cache")

    assert first.cache_dir == second.cache_dir
    assert called["count"] == 0
    assert second.assembly is not None and second.assembly.exists()


def test_parse_vcf_records_ignores_headers(tmp_path: Path):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "ref\t2\t.\tA\tG\t.\tPASS\t.\n"
    )
    assert parse_vcf_records(vcf) == [VariantRecord("ref", 2, "A", "G")]
