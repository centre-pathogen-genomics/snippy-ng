from pathlib import Path

import pytest

import snippy_ng.pipelines as _pl
from snippy_ng.stages.masks import DepthBedsFromBam, DepthMaskFromBed
from snippy_ng.stages.vcf import AddDeletionstoVCF
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
        "--outdir", str(paths["out"]),
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert DepthBedsFromBam in _stage_types()

    depth_beds = _stage_instances(DepthBedsFromBam)[0]
    add_deletions = _stage_instances(AddDeletionstoVCF)[0]
    depth_mask = _stage_instances(DepthMaskFromBed)[0]

    assert add_deletions.zero_depth_bed == Path("snippy.zerodepth.bed")
    assert add_deletions.zero_depth_bed == depth_beds.output.zero_depth_bed
    assert depth_mask.mask_bed == Path("snippy.mindepth.bed")
    assert depth_mask.mask_bed == depth_beds.output.min_depth_bed


def test_short_pipeline_depth_mask_zero_keeps_zero_depth_deletions_only(tmp_path):
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
        "--depth-mask", "0",
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert DepthBedsFromBam in _stage_types()
    assert AddDeletionstoVCF in _stage_types()
    assert DepthMaskFromBed not in _stage_types()


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
        "--outdir", str(paths["out"]),
        "--caller", "freebayes",
        "--skip-check",
    ])

    assert_cli_result(result, 0, True)
    assert DepthBedsFromBam in _stage_types()

    depth_beds = _stage_instances(DepthBedsFromBam)[0]
    add_deletions = _stage_instances(AddDeletionstoVCF)[0]
    depth_mask = _stage_instances(DepthMaskFromBed)[0]

    assert add_deletions.zero_depth_bed == depth_beds.output.zero_depth_bed
    assert depth_mask.mask_bed == depth_beds.output.min_depth_bed
