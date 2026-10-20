from pathlib import Path

import pytest

from snippy_ng.context import Context
from snippy_ng.exceptions import MSAValidationError
from snippy_ng.pipelines.core import CorePipelineBuilder
from snippy_ng.stages.alignment_filter import (
    CheckAlignmentClustersByPipelineType,
    FilterAlignmentByAlignedPercentage,
)
from snippy_ng.stages.core import (
    CombineFastaFile,
    DistleDistanceMatrix,
    SoftCoreError,
    SoftCoreFilter,
)


def _build_core_pipeline(tmp_path, mask=None):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "reference.fa").write_text(">ref\nAAAA\n")
    (reference / "reference.fa.fai").write_text("ref\t4\t0\t4\t5\n")
    (reference / "reference.dict").write_text("ref\t4\n")
    (reference / "reference.gff").write_text("##gff-version 3\n")
    (reference / "metadata.json").write_text('{"prefix": "reference"}')

    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()
    (snippy_dir / "snippy.qc.tsv").write_text("sample\tpipeline_type\nsample1\tshort\n")

    return CorePipelineBuilder(
        snippy_dirs=[snippy_dir],
        reference=reference,
        prefix="core",
        mask=mask,
    ).build()


def test_core_pipeline_soft_core_uses_filtered_alignment(tmp_path):
    pipeline = _build_core_pipeline(tmp_path)

    soft_core_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, SoftCoreFilter)
    )
    filter_stage = next(
        stage
        for stage in pipeline.stages
        if isinstance(stage, FilterAlignmentByAlignedPercentage)
    )
    technical_check = next(
        stage
        for stage in pipeline.stages
        if isinstance(stage, CheckAlignmentClustersByPipelineType)
    )

    assert soft_core_stage.aln == filter_stage.output.filtered_aln
    assert filter_stage.output.filtered_aln in pipeline.outputs_to_keep
    assert filter_stage.output.filter_stats in pipeline.outputs_to_keep
    assert filter_stage.remove_outliers is False
    assert filter_stage.outlier_threshold == 0.20
    assert technical_check.filter_stats == filter_stage.output.filter_stats
    assert technical_check.qc_files == [tmp_path / "sample1" / "snippy.qc.tsv"]
    assert technical_check.output.summary_tsv in pipeline.outputs_to_keep


def test_core_pipeline_reference_id_uses_metadata_prefix_when_distinct(tmp_path):
    pipeline = _build_core_pipeline(tmp_path)

    combine_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile)
    )

    assert combine_stage.reference_id == "reference"


def test_core_pipeline_reference_id_avoids_sample_name_collision(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "reference.fa").write_text(">ref\nAAAA\n")
    (reference / "reference.fa.fai").write_text("ref\t4\t0\t4\t5\n")
    (reference / "reference.dict").write_text("ref\t4\n")
    (reference / "reference.gff").write_text("##gff-version 3\n")
    (reference / "metadata.json").write_text('{"prefix": "AUSMDU00112379"}')

    snippy_dir = tmp_path / "AUSMDU00112379"
    snippy_dir.mkdir()
    (snippy_dir / "snippy.qc.tsv").write_text(
        "sample\tpipeline_type\nAUSMDU00112379\tshort\n"
    )

    pipeline = CorePipelineBuilder(
        snippy_dirs=[snippy_dir],
        reference=reference,
        prefix="core",
    ).build()

    combine_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile)
    )

    assert combine_stage.reference_id == "reference"


def test_core_pipeline_reference_id_suffixes_if_reference_is_also_a_sample(tmp_path):
    reference = tmp_path / "prepared_reference"
    reference.mkdir()
    (reference / "reference.fa").write_text(">ref\nAAAA\n")
    (reference / "reference.fa.fai").write_text("ref\t4\t0\t4\t5\n")
    (reference / "reference.dict").write_text("ref\t4\n")
    (reference / "reference.gff").write_text("##gff-version 3\n")
    (reference / "metadata.json").write_text('{"prefix": "AUSMDU00112379"}')

    snippy_dirs = []
    for sample_name in ("AUSMDU00112379", "reference"):
        snippy_dir = tmp_path / sample_name
        snippy_dir.mkdir()
        (snippy_dir / "snippy.qc.tsv").write_text(
            f"sample\tpipeline_type\n{sample_name}\tshort\n"
        )
        snippy_dirs.append(snippy_dir)

    pipeline = CorePipelineBuilder(
        snippy_dirs=snippy_dirs,
        reference=reference,
        prefix="core",
    ).build()

    combine_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile)
    )

    assert combine_stage.reference_id == "reference_2"


def _combine_with_mask(tmp_path, bed_text, mask_char="N"):
    reference = tmp_path / "ref.fa"
    reference.write_text(">c1\nAAAAAA\n>c2\nCCCCCC\n")
    sample_dir = tmp_path / "s1"
    sample_dir.mkdir()
    (sample_dir / "s1.fna").write_text(">c1\nACGTAC\n>c2\nGGTTGG\n")
    bed = tmp_path / "mask.bed"
    bed.write_text(bed_text)

    stage = CombineFastaFile(
        snippy_dirs=[sample_dir],
        reference=reference,
        mask_bed=bed,
        mask_char=mask_char,
        prefix=str(tmp_path / "core"),
    )
    stage.build_concatenated_alignment([sample_dir], reference, stage.output.aln)
    records = {}
    name = None
    for line in stage.output.aln.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
        else:
            records[name] = line
    return records


def test_combine_fasta_file_masks_reference_and_samples_using_contig_offsets(tmp_path):
    records = _combine_with_mask(
        tmp_path,
        "# comment\n"
        "track name=mask\n"
        "c1\t1\t3\n"
        "c2\t0\t2\textra\tcolumns\n"
        "c2\t1\t3\n",  # overlaps previous interval and must be merged
    )

    assert records["reference"] == "ANNAAA" + "NNNCCC"
    assert records["s1"] == "ANNTAC" + "NNNTGG"


def test_combine_fasta_file_mask_char_is_configurable(tmp_path):
    records = _combine_with_mask(tmp_path, "c1\t0\t1\n", mask_char="-")

    assert records["s1"] == "-CGTACGGTTGG"
    assert records["reference"] == "-AAAAACCCCCC"


@pytest.mark.parametrize(
    "bed_text, message",
    [
        ("c3\t0\t1\n", "not in the reference"),
        ("c1\t0\t7\n", "invalid for contig"),
        ("c1\t4\t2\n", "invalid for contig"),
        ("c1\tfoo\t2\n", "must be integers"),
        ("c1\t1\n", "at least 3 BED columns"),
    ],
)
def test_combine_fasta_file_rejects_invalid_mask(tmp_path, bed_text, message):
    with pytest.raises(MSAValidationError, match=message):
        _combine_with_mask(tmp_path, bed_text)


def test_core_pipeline_passes_mask_to_combine_stage(tmp_path):
    mask = tmp_path / "mask.bed"
    mask.write_text("ref\t0\t1\n")
    pipeline = _build_core_pipeline(tmp_path, mask=mask)

    combine_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile)
    )

    assert combine_stage.mask_bed == mask


def test_core_pipeline_has_no_mask_by_default(tmp_path):
    pipeline = _build_core_pipeline(tmp_path)

    combine_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile)
    )

    assert combine_stage.mask_bed is None


def test_distle_distance_matrix_uses_tabular_output(tmp_path):
    stage = DistleDistanceMatrix(
        aln=tmp_path / "core.full.aln",
        prefix="core.full",
    )

    commands = stage.create_commands(Context())

    assert len(commands) == 1
    assert commands[0].command == [
        "distle",
        "--threads",
        str(Context().cpus),
        "-o",
        "tabular",
        str(tmp_path / "core.full.aln"),
        str(tmp_path / "core.full.distance.tsv"),
    ]
    assert stage.output.phylip == tmp_path / "core.full.distance.tsv"


def test_core_pipeline_adds_distle_for_full_and_soft_core_alignments(tmp_path):
    pipeline = _build_core_pipeline(tmp_path)

    distle_stages = [
        stage for stage in pipeline.stages if isinstance(stage, DistleDistanceMatrix)
    ]
    combine_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile)
    )
    soft_core_stage = next(
        stage for stage in pipeline.stages if isinstance(stage, SoftCoreFilter)
    )

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
