import pytest
from snippy_ng.stages.reporting import SampleReport
from snippy_ng.stages.stats import FastaCompositionStats, SampleQcSummary
from snippy_ng.stages.vcf import VcfPassFilter, VcfToTab
from tests.cli.helpers import apply_cli_case_overrides, assert_cli_result, get_bad_reference_target, run_cli_command, write_dummy_files


@pytest.fixture(autouse=True)
def stub_everything(stub_pipeline, stub_reference_format, stub_common_stages, stub_asm_stages):
    """
    Combine all required stubs for assembly CLI tests.
    This fixture runs automatically for every test.
    """
    pass


##############################################################################
#                               2.  TEST DATA                                 #
##############################################################################

@pytest.mark.parametrize(
    "case_name, extra, expect_exit, expect_run",
    [
        (
            "assembly_ok",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "assembly_ok_nucmer",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--caller",   "nucmer",
                "--outdir",    p["out"],
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "assembly_ok_minimap_asm5",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--caller",   "paftools",
                "--outdir",    p["out"],
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "check_only",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--check",
                "--skip-check",
            ],
            0,
            False,
        ),
        (
            "outdir_exists",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            2,
            False,
        ),
        (
            "bad_reference",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            1,
            False,
        ),
    ],
)
def test_asm_cli(monkeypatch, tmp_path, case_name, extra, expect_exit, expect_run):
    """
    Parameterised test for the `asm` command.
    """

    # --------------- Arrange --------------------------------------------------
    paths = {
        "ref": tmp_path / "ref.fa",
        "asm": tmp_path / "assembly.fa",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "asm"])
    apply_cli_case_overrides(
        monkeypatch,
        case_name,
        paths["out"],
        bad_reference_target=get_bad_reference_target("asm"),
    )

    args = ["asm"] + extra(paths)
    result = run_cli_command(args)
    assert_cli_result(result, expect_exit, expect_run)


def test_asm_pipeline_includes_vcf_only_sample_report(monkeypatch, tmp_path):
    from snippy_ng.pipelines.asm import AsmPipelineBuilder

    paths = {
        "ref": tmp_path / "ref.fa",
        "asm": tmp_path / "assembly.fa",
    }
    write_dummy_files(paths, ["ref", "asm"])

    pipeline = AsmPipelineBuilder(
        reference=paths["ref"],
        assembly=paths["asm"],
        prefix="sample",
    ).build()

    sample_report_stages = [stage for stage in pipeline.stages if isinstance(stage, SampleReport)]
    pass_filter = next(stage for stage in pipeline.stages if isinstance(stage, VcfPassFilter))
    variants_tab = next(stage for stage in pipeline.stages if isinstance(stage, VcfToTab))
    fasta_qc = next(stage for stage in pipeline.stages if isinstance(stage, FastaCompositionStats))
    sample_qc = next(stage for stage in pipeline.stages if isinstance(stage, SampleQcSummary))
    assert len(sample_report_stages) == 1
    assert sample_report_stages[0].alignment is None
    assert sample_report_stages[0].reference is None
    assert sample_report_stages[0].variant_scope == "all"
    assert variants_tab.vcf == pass_filter.output.vcf
    assert variants_tab.output.tab in pipeline.outputs_to_keep
    assert fasta_qc.output.summary_tsv not in pipeline.outputs_to_keep
    assert sample_qc.pipeline_type == "asm"
    assert sample_qc.reads_tsv is None
    assert sample_qc.alignment_tsv is None
    assert sample_qc.output.qc_tsv in pipeline.outputs_to_keep


def test_asm_pipeline_can_omit_sample_report(monkeypatch, tmp_path):
    from snippy_ng.pipelines.asm import AsmPipelineBuilder

    paths = {
        "ref": tmp_path / "ref.fa",
        "asm": tmp_path / "assembly.fa",
    }
    write_dummy_files(paths, ["ref", "asm"])

    pipeline = AsmPipelineBuilder(
        reference=paths["ref"],
        assembly=paths["asm"],
        report=False,
        prefix="sample",
    ).build()

    assert not any(isinstance(stage, SampleReport) for stage in pipeline.stages)


def test_asm_pipeline_adds_context_filter_when_enabled(monkeypatch, tmp_path):
    from snippy_ng.pipelines.asm import AsmPipelineBuilder
    from snippy_ng.stages.vcf import VariantContextFilter

    monkeypatch.setenv("SNIPPY_NG_VARIANT_CONTEXT_MIN_SNP_DISTANCE_TO_INDEL", "10")
    paths = {
        "ref": tmp_path / "ref.fa",
        "asm": tmp_path / "assembly.fa",
    }
    write_dummy_files(paths, ["ref", "asm"])

    pipeline = AsmPipelineBuilder(
        reference=paths["ref"],
        assembly=paths["asm"],
        prefix="sample",
    ).build()

    assert any(isinstance(stage, VariantContextFilter) for stage in pipeline.stages)
