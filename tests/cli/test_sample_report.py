from pathlib import Path

import snippy_ng.pipelines.reports as report_pipeline_module
from snippy_ng.stages.reporting import SampleReport
from tests.cli.conftest import DummyPipeline
from tests.cli.helpers import make_prepared_reference, run_cli_command, stub_load_or_prepare_reference


def test_sample_report_cli_builds_pipeline(monkeypatch, tmp_path):
    monkeypatch.setattr(
        report_pipeline_module.pipelines,
        "SnippyPipeline",
        lambda stages=None, outputs_to_keep=None: DummyPipeline(stages=stages, outputs_to_keep=outputs_to_keep),
    )

    vcf = tmp_path / "sample.vcf"
    alignment = tmp_path / "sample.cram"
    reference = tmp_path / "ref.fa"
    reference_index = tmp_path / "ref.fa.fai"
    outdir = tmp_path / "report"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    alignment.write_text("cram")
    reference.write_text(">chr1\nA\n")
    reference_index.write_text("chr1\t1\t6\t1\t2\n")

    result = run_cli_command(
        [
            "utils",
            "sample-report",
            str(vcf),
            "--alignment",
            str(alignment),
            "--reference",
            str(reference),
            "--variant-scope",
            "all",
            "--window-size",
            "500",
            "--outdir",
            str(outdir),
            "--skip-check",
        ]
    )

    assert result.exit_code == 0, result.output
    pipeline = DummyPipeline.last
    assert pipeline.ran is True
    stage = next(stage for stage in pipeline.stages if isinstance(stage, SampleReport))
    assert stage.variant_scope == "all"
    assert stage.window_size == 500
    assert stage.reference_index == Path("reference/ref.fa.fai")


def test_sample_report_cli_builds_vcf_only_pipeline(monkeypatch, tmp_path):
    monkeypatch.setattr(
        report_pipeline_module.pipelines,
        "SnippyPipeline",
        lambda stages=None, outputs_to_keep=None: DummyPipeline(stages=stages, outputs_to_keep=outputs_to_keep),
    )

    vcf = tmp_path / "sample.vcf"
    outdir = tmp_path / "report"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    result = run_cli_command(
        [
            "utils",
            "sample-report",
            str(vcf),
            "--outdir",
            str(outdir),
            "--skip-check",
        ]
    )

    assert result.exit_code == 0, result.output
    pipeline = DummyPipeline.last
    assert pipeline.ran is True
    assert len(pipeline.stages) == 1
    stage = pipeline.stages[0]
    assert isinstance(stage, SampleReport)
    assert stage.alignment is None
    assert stage.reference is None


def test_sample_report_cli_requires_reference_with_alignment(tmp_path):
    vcf = tmp_path / "sample.vcf"
    alignment = tmp_path / "sample.cram"
    vcf.write_text("##fileformat=VCFv4.2\n")
    alignment.write_text("cram")

    result = run_cli_command(
        [
            "utils",
            "sample-report",
            str(vcf),
            "--alignment",
            str(alignment),
            "--outdir",
            str(tmp_path / "report"),
            "--skip-check",
        ]
    )

    assert result.exit_code == 2
    assert "--reference is required when --alignment is provided" in result.output


def test_short_pipeline_includes_sample_report_by_default(monkeypatch, tmp_path):
    from snippy_ng.pipelines.short import ShortPipelineBuilder

    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.short.load_or_prepare_reference",
    )
    bam = tmp_path / "reads.bam"
    bam.write_text("bam")

    pipeline = ShortPipelineBuilder(
        reference=ref_file,
        reads=[],
        bam=bam,
        prefix="sample",
    ).build()

    sample_report_stages = [stage for stage in pipeline.stages if isinstance(stage, SampleReport)]
    assert len(sample_report_stages) == 1
    assert sample_report_stages[0].variant_scope == "pass"


def test_long_pipeline_can_omit_sample_report(monkeypatch, tmp_path):
    from snippy_ng.pipelines.long import LongPipelineBuilder

    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.long.load_or_prepare_reference",
    )
    bam = tmp_path / "reads.bam"
    bam.write_text("bam")

    pipeline = LongPipelineBuilder(
        reference=ref_file,
        reads=None,
        bam=bam,
        caller="freebayes",
        sample_report=False,
        prefix="sample",
    ).build()

    assert not any(isinstance(stage, SampleReport) for stage in pipeline.stages)
