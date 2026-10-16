from pathlib import Path

import pytest

import snippy_ng.pipelines as _pl
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.filtering import BamReferenceValidator, SamtoolsFilter
from snippy_ng.stages.masks import DepthBedsFromBam, ApplyDepthMaskToFasta, MinDepthThresholdFromBam
from snippy_ng.stages.vcf import AddDeletionsToVCF, VcfFilterLong, VcfFilterShort
from tests.cli.helpers import assert_cli_result, run_cli_command, write_dummy_files


@pytest.fixture(autouse=True)
def stub_everything(stub_pipeline, stub_reference_format, stub_common_stages, stub_short_stages, stub_long_stages):
    pass


def _stage_types():
    return [type(stage) for stage in _pl.SnippyPipeline.last.stages]


def _stage_instances(stage_type):
    return [stage for stage in _pl.SnippyPipeline.last.stages if isinstance(stage, stage_type)]


def test_short_pipeline_uses_combined_depth_beds(tmp_path):
    paths = {
        "ref": tmp_path / "ref.fa",
        "r1": tmp_path / "reads_1.fq",
        "r2": tmp_path / "reads_2.fq",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "r1", "r2"])

    result = run_cli_command([
        "short",
        "--reference", str(paths["ref"]),
        "--R1", str(paths["r1"]),
        "--R2", str(paths["r2"]),
        "--min-depth", "10",
        "--outdir", str(paths["out"]),
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert DepthBedsFromBam in _stage_types()
    assert MinDepthThresholdFromBam not in _stage_types()

    depth_beds = _stage_instances(DepthBedsFromBam)[0]
    add_deletions = _stage_instances(AddDeletionsToVCF)[0]
    depth_mask = _stage_instances(ApplyDepthMaskToFasta)[0]
    vcf_filter = _stage_instances(VcfFilterShort)[0]

    assert add_deletions.zero_depth_bed == Path("snippy.zerodepth.bed")
    assert add_deletions.zero_depth_bed == depth_beds.output.zero_depth_bed
    assert depth_mask.mask_bed == Path("snippy.mindepth.bed")
    assert depth_mask.mask_bed == depth_beds.output.min_depth_bed
    # fixed --min-depth is used directly; no dynamic threshold stage
    assert vcf_filter.min_depth == 10
    assert vcf_filter.min_depth_file is None
    assert depth_beds.min_depth == 10
    assert depth_beds.min_depth_file is None
    assert depth_mask.min_depth_file is None
    assert _stage_types().index(BcftoolsPseudoAlignment) < _stage_types().index(ApplyDepthMaskToFasta)


def test_short_pipeline_min_depth_zero_keeps_zero_depth_deletions_only(tmp_path):
    paths = {
        "ref": tmp_path / "ref.fa",
        "r1": tmp_path / "reads_1.fq",
        "r2": tmp_path / "reads_2.fq",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "r1", "r2"])

    result = run_cli_command([
        "short",
        "--reference", str(paths["ref"]),
        "--R1", str(paths["r1"]),
        "--R2", str(paths["r2"]),
        "--outdir", str(paths["out"]),
        "--min-depth", "0",
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert DepthBedsFromBam in _stage_types()
    assert AddDeletionsToVCF in _stage_types()
    assert ApplyDepthMaskToFasta not in _stage_types()
    assert MinDepthThresholdFromBam not in _stage_types()


def test_short_pipeline_validates_user_bam_before_filtering(tmp_path):
    paths = {
        "ref": tmp_path / "ref.fa",
        "bam": tmp_path / "reads.bam",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "bam"])

    result = run_cli_command([
        "short",
        "--reference", str(paths["ref"]),
        "--bam", str(paths["bam"]),
        "--outdir", str(paths["out"]),
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert _stage_types().index(BamReferenceValidator) < _stage_types().index(SamtoolsFilter)


def test_long_pipeline_uses_combined_depth_beds(tmp_path):
    paths = {
        "ref": tmp_path / "ref.fa",
        "reads": tmp_path / "reads.fq",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "reads"])

    result = run_cli_command([
        "long",
        "--reference", str(paths["ref"]),
        "--reads", str(paths["reads"]),
        "--min-depth", "10",
        "--outdir", str(paths["out"]),
        "--caller", "freebayes",
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert DepthBedsFromBam in _stage_types()
    assert MinDepthThresholdFromBam not in _stage_types()

    depth_beds = _stage_instances(DepthBedsFromBam)[0]
    add_deletions = _stage_instances(AddDeletionsToVCF)[0]
    depth_mask = _stage_instances(ApplyDepthMaskToFasta)[0]
    vcf_filter = _stage_instances(VcfFilterLong)[0]

    assert add_deletions.zero_depth_bed == depth_beds.output.zero_depth_bed
    assert depth_mask.mask_bed == depth_beds.output.min_depth_bed
    # fixed --min-depth is used directly; no dynamic threshold stage
    assert vcf_filter.min_depth == 10
    assert vcf_filter.min_depth_file is None
    assert depth_beds.min_depth_file is None
    assert depth_mask.min_depth_file is None
    assert _stage_types().index(BcftoolsPseudoAlignment) < _stage_types().index(ApplyDepthMaskToFasta)


def test_long_pipeline_validates_user_bam_before_filtering(tmp_path):
    paths = {
        "ref": tmp_path / "ref.fa",
        "bam": tmp_path / "reads.bam",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "bam"])

    result = run_cli_command([
        "long",
        "--reference", str(paths["ref"]),
        "--bam", str(paths["bam"]),
        "--outdir", str(paths["out"]),
        "--caller", "freebayes",
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert _stage_types().index(BamReferenceValidator) < _stage_types().index(SamtoolsFilter)


def test_short_pipeline_default_min_depth_is_dynamic(tmp_path):
    paths = {
        "ref": tmp_path / "ref.fa",
        "r1": tmp_path / "reads_1.fq",
        "r2": tmp_path / "reads_2.fq",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "r1", "r2"])

    result = run_cli_command([
        "short",
        "--reference", str(paths["ref"]),
        "--R1", str(paths["r1"]),
        "--R2", str(paths["r2"]),
        "--outdir", str(paths["out"]),
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert MinDepthThresholdFromBam in _stage_types()

    threshold_txt = _stage_instances(MinDepthThresholdFromBam)[0].output.threshold_txt
    depth_beds = _stage_instances(DepthBedsFromBam)[0]
    vcf_filter = _stage_instances(VcfFilterShort)[0]
    depth_mask = _stage_instances(ApplyDepthMaskToFasta)[0]

    assert vcf_filter.min_depth_file == threshold_txt
    assert depth_beds.min_depth_file == threshold_txt
    assert depth_mask.min_depth_file == threshold_txt
    # threshold must be computed before any stage consumes it
    threshold_index = _stage_types().index(MinDepthThresholdFromBam)
    assert threshold_index < _stage_types().index(VcfFilterShort)
    assert threshold_index < _stage_types().index(DepthBedsFromBam)


def test_long_pipeline_default_min_depth_is_dynamic(tmp_path):
    paths = {
        "ref": tmp_path / "ref.fa",
        "reads": tmp_path / "reads.fq",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "reads"])

    result = run_cli_command([
        "long",
        "--reference", str(paths["ref"]),
        "--reads", str(paths["reads"]),
        "--outdir", str(paths["out"]),
        "--caller", "freebayes",
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert MinDepthThresholdFromBam in _stage_types()

    threshold_txt = _stage_instances(MinDepthThresholdFromBam)[0].output.threshold_txt
    depth_beds = _stage_instances(DepthBedsFromBam)[0]
    vcf_filter = _stage_instances(VcfFilterLong)[0]
    depth_mask = _stage_instances(ApplyDepthMaskToFasta)[0]

    assert vcf_filter.min_depth_file == threshold_txt
    assert depth_beds.min_depth_file == threshold_txt
    assert depth_mask.min_depth_file == threshold_txt
    # threshold must be computed before any stage consumes it
    threshold_index = _stage_types().index(MinDepthThresholdFromBam)
    assert threshold_index < _stage_types().index(VcfFilterLong)
    assert threshold_index < _stage_types().index(DepthBedsFromBam)
