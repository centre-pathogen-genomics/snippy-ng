from pathlib import Path

import pytest

from snippy_ng.context import Context
from snippy_ng.pipelines.core import CorePipelineBuilder
from snippy_ng.stages.alignment_filter import FilterAlignmentByAlignedPercentage
from snippy_ng.stages.core import CombineFastaFile, DistleDistanceMatrix, SoftCoreError, SoftCoreFilter


def _build_core_pipeline(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "reference.fa").write_text(">ref\nAAAA\n")
    (reference / "reference.fa.fai").write_text("ref\t4\t0\t4\t5\n")
    (reference / "reference.dict").write_text("ref\t4\n")
    (reference / "reference.gff").write_text("##gff-version 3\n")
    (reference / "metadata.json").write_text('{"prefix": "reference"}')

    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()

    return CorePipelineBuilder(
        snippy_dirs=[snippy_dir],
        reference=reference,
        prefix="core",
    ).build()


def test_core_pipeline_soft_core_uses_filtered_alignment(tmp_path):
    pipeline = _build_core_pipeline(tmp_path)

    soft_core_stage = next(stage for stage in pipeline.stages if isinstance(stage, SoftCoreFilter))
    filter_stage = next(stage for stage in pipeline.stages if isinstance(stage, FilterAlignmentByAlignedPercentage))

    assert soft_core_stage.aln == filter_stage.output.filtered_aln
    assert filter_stage.output.filtered_aln in pipeline.outputs_to_keep
    assert filter_stage.output.filter_stats in pipeline.outputs_to_keep
    assert filter_stage.inclusion_threshold == 0.20


def test_distle_distance_matrix_uses_tabular_output(tmp_path):
    stage = DistleDistanceMatrix(
        aln=tmp_path / "core.full.aln",
        prefix="core.full",
    )

    commands = stage.create_commands(Context())

    assert len(commands) == 1
    assert commands[0].command == [
        "distle",
        "--threads", str(Context().cpus),
        "-o",
        "tabular",
        str(tmp_path / "core.full.aln"),
        str(tmp_path / "core.full.distance.tsv"),
    ]
    assert stage.output.phylip == tmp_path / "core.full.distance.tsv"


def test_core_pipeline_adds_distle_for_full_and_soft_core_alignments(tmp_path):
    pipeline = _build_core_pipeline(tmp_path)

    distle_stages = [stage for stage in pipeline.stages if isinstance(stage, DistleDistanceMatrix)]
    combine_stage = next(stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile))
    soft_core_stage = next(stage for stage in pipeline.stages if isinstance(stage, SoftCoreFilter))

    assert len(distle_stages) == 2
    assert distle_stages[0].aln == combine_stage.output.aln
    assert distle_stages[0].output.phylip == Path("core.full.distance.tsv")
    assert distle_stages[1].aln == soft_core_stage.output.soft_core
    assert distle_stages[1].output.phylip == Path("core.095.distance.tsv")


def test_soft_core_empty_error_describes_missing_aligned_calls(tmp_path):
    stage = SoftCoreFilter(
        aln=tmp_path / "core.filtered.aln",
        core_threshold=0.95,
        prefix=str(tmp_path / "core"),
    )
    stage.output.soft_core.write_text(">reference\n\n")

    with pytest.raises(
        SoftCoreError,
        match="too few aligned or called bases",
    ):
        stage.test_soft_core_is_not_empty()
