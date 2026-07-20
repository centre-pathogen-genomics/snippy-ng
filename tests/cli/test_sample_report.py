from pathlib import Path

import snippy_ng.pipelines.reports as report_pipeline_module
from snippy_ng.stages.copy import CopyFile
from snippy_ng.stages.downsample import SamtoolsDownsampleAlignment
from snippy_ng.stages.reporting import SampleReport
from snippy_ng.stages.stats import FastaCompositionStats, SampleQcSummary, SamtoolsAlignmentQcStats, VcfStats
from snippy_ng.stages.vcf import CollapseDiploidGenotypes, VcfPassFilter, VcfToTab
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
            "report",
            "sample",
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
    assert stage.exclude_supplementary is False
    assert not any(isinstance(stage, SamtoolsDownsampleAlignment) for stage in pipeline.stages)


def test_sample_report_cli_can_downsample_alignment(monkeypatch, tmp_path):
    monkeypatch.setattr(
        report_pipeline_module.pipelines,
        "SnippyPipeline",
        lambda stages=None, outputs_to_keep=None: DummyPipeline(stages=stages, outputs_to_keep=outputs_to_keep),
    )

    vcf = tmp_path / "sample.vcf"
    alignment = tmp_path / "sample.bam"
    reference = tmp_path / "ref.fa"
    outdir = tmp_path / "report"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    alignment.write_text("bam")
    reference.write_text(">chr1\nA\n")
    stub_load_or_prepare_reference(
        monkeypatch,
        reference,
        target="snippy_ng.pipelines.reports.load_or_prepare_reference",
    )

    result = run_cli_command(
        [
            "utils",
            "report",
            "sample",
            str(vcf),
            "--alignment",
            str(alignment),
            "--reference",
            str(reference),
            "--downsample",
            "0.25",
            "--outdir",
            str(outdir),
            "--skip-check",
        ]
    )

    assert result.exit_code == 0, result.output
    pipeline = DummyPipeline.last
    downsample_stage = next(stage for stage in pipeline.stages if isinstance(stage, SamtoolsDownsampleAlignment))
    report_stage = next(stage for stage in pipeline.stages if isinstance(stage, SampleReport))
    assert downsample_stage.alignment == alignment
    assert downsample_stage.fraction == 0.25
    assert report_stage.alignment == downsample_stage.output.bam


def test_sample_report_cli_rejects_downsample_without_alignment(tmp_path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")

    result = run_cli_command(
        [
            "utils",
            "report",
            "sample",
            str(vcf),
            "--downsample",
            "0.25",
            "--outdir",
            str(tmp_path / "report"),
            "--skip-check",
        ]
    )

    assert result.exit_code == 2
    assert "--alignment is required when --downsample is provided" in result.output


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
            "report",
            "sample",
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
            "report",
            "sample",
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
    assert sample_report_stages[0].variant_scope == "all"


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
        report=False,
        prefix="sample",
    ).build()

    assert not any(isinstance(stage, SampleReport) for stage in pipeline.stages)


def test_short_pipeline_uses_final_vcf_copy_downstream(monkeypatch, tmp_path):
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
        report_scope="all",
        prefix="sample",
    ).build()

    collapse_stage = next(stage for stage in pipeline.stages if isinstance(stage, CollapseDiploidGenotypes))
    final_vcf = next(stage for stage in pipeline.stages if isinstance(stage, CopyFile) and stage.output_path == Path("sample.all.vcf"))
    vcf_stats = next(stage for stage in pipeline.stages if isinstance(stage, VcfStats))
    pass_filter = next(stage for stage in pipeline.stages if isinstance(stage, VcfPassFilter))
    variants_tab = next(stage for stage in pipeline.stages if isinstance(stage, VcfToTab))
    sample_report = next(stage for stage in pipeline.stages if isinstance(stage, SampleReport))
    alignment_qc = next(stage for stage in pipeline.stages if isinstance(stage, SamtoolsAlignmentQcStats))
    fasta_qc = next(stage for stage in pipeline.stages if isinstance(stage, FastaCompositionStats))
    sample_qc = next(stage for stage in pipeline.stages if isinstance(stage, SampleQcSummary))

    assert final_vcf.input == collapse_stage.output.vcf
    assert vcf_stats.vcf == final_vcf.output.copied_file
    assert pass_filter.vcf == final_vcf.output.copied_file
    assert variants_tab.vcf == pass_filter.output.vcf
    assert variants_tab.output.tab in pipeline.outputs_to_keep
    assert sample_report.vcf == final_vcf.output.copied_file
    assert sample_report.exclude_supplementary is False
    assert alignment_qc.bam == Path("sample.filtered.bam")
    assert fasta_qc.fasta == Path("sample.fna")
    assert sample_qc.output.qc_tsv in pipeline.outputs_to_keep
    assert alignment_qc.output.summary_tsv not in pipeline.outputs_to_keep
    assert fasta_qc.output.summary_tsv not in pipeline.outputs_to_keep


def test_long_clair3_report_excludes_supplementary_alignments(monkeypatch, tmp_path):
    from snippy_ng.pipelines.long import LongPipelineBuilder

    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.long.load_or_prepare_reference",
    )
    bam = tmp_path / "reads.bam"
    bam.write_text("bam")
    clair3_model = tmp_path / "clair3-model"
    clair3_model.mkdir()

    pipeline = LongPipelineBuilder(
        reference=ref_file,
        reads=None,
        bam=bam,
        caller="clair3",
        model=clair3_model,
        prefix="sample",
    ).build()

    sample_report = next(stage for stage in pipeline.stages if isinstance(stage, SampleReport))

    assert sample_report.exclude_supplementary is True


def test_long_pipeline_uses_final_vcf_copy_downstream(monkeypatch, tmp_path):
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
        report_scope="all",
        prefix="sample",
    ).build()

    collapse_stage = next(stage for stage in pipeline.stages if isinstance(stage, CollapseDiploidGenotypes))
    final_vcf = next(stage for stage in pipeline.stages if isinstance(stage, CopyFile) and stage.output_path == Path("sample.all.vcf"))
    vcf_stats = next(stage for stage in pipeline.stages if isinstance(stage, VcfStats))
    pass_filter = next(stage for stage in pipeline.stages if isinstance(stage, VcfPassFilter))
    variants_tab = next(stage for stage in pipeline.stages if isinstance(stage, VcfToTab))
    sample_report = next(stage for stage in pipeline.stages if isinstance(stage, SampleReport))
    alignment_qc = next(stage for stage in pipeline.stages if isinstance(stage, SamtoolsAlignmentQcStats))
    fasta_qc = next(stage for stage in pipeline.stages if isinstance(stage, FastaCompositionStats))
    sample_qc = next(stage for stage in pipeline.stages if isinstance(stage, SampleQcSummary))

    assert final_vcf.input == collapse_stage.output.vcf
    assert vcf_stats.vcf == final_vcf.output.copied_file
    assert pass_filter.vcf == final_vcf.output.copied_file
    assert variants_tab.vcf == pass_filter.output.vcf
    assert variants_tab.output.tab in pipeline.outputs_to_keep
    assert sample_report.vcf == final_vcf.output.copied_file
    assert alignment_qc.bam == Path("sample.filtered.bam")
    assert fasta_qc.fasta == Path("sample.fna")
    assert sample_qc.output.qc_tsv in pipeline.outputs_to_keep
    assert alignment_qc.output.summary_tsv not in pipeline.outputs_to_keep
    assert fasta_qc.output.summary_tsv not in pipeline.outputs_to_keep


def test_short_pipeline_skips_genotype_collapse_when_haploid_disabled(monkeypatch, tmp_path):
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
        haploid=False,
        prefix="sample",
    ).build()

    assert not any(isinstance(stage, CollapseDiploidGenotypes) for stage in pipeline.stages)
    final_vcf = next(stage for stage in pipeline.stages if isinstance(stage, CopyFile) and stage.output_path == Path("sample.all.vcf"))
    assert final_vcf.input == Path("sample.annotated.vcf")


def test_long_pipeline_skips_genotype_collapse_when_haploid_disabled(monkeypatch, tmp_path):
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
        haploid=False,
        prefix="sample",
    ).build()

    assert not any(isinstance(stage, CollapseDiploidGenotypes) for stage in pipeline.stages)
    final_vcf = next(stage for stage in pipeline.stages if isinstance(stage, CopyFile) and stage.output_path == Path("sample.all.vcf"))
    assert final_vcf.input == Path("sample.annotated.vcf")
