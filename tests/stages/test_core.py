import csv
from pathlib import Path

from Bio import SeqIO

from snippy_ng.logging import logger
from snippy_ng.pipelines.core import CorePipelineBuilder
from snippy_ng.stages.core import CombineFastaFile, DistleDistanceMatrix, FilterAlignmentByAlignedPercentage, SoftCoreFilter
from snippy_ng.context import Context


def test_filter_alignment_by_aligned_percentage_filters_outlier_sample(tmp_path):
    aln = tmp_path / "core.full.aln"
    aln.write_text(
        ">reference\nAAAAAA\n"
        ">sample_a\nAAAAAA\n"
        ">sample_b\nAAA-AAA\n"
        ">sample_c\nAA----\n"
    )
    aligned_tsv = tmp_path / "core.aligned.tsv"
    aligned_tsv.write_text(
        "sequence\taligned\n"
        "reference\t100.00\n"
        "sample_a\t100.00\n"
        "sample_b\t85.00\n"
        "sample_c\t20.00\n"
    )
    filtered_aln = tmp_path / "core.filtered.aln"

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference", "sample_a", "sample_b"]

    with aligned_tsv.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0] == {
        "sequence": "reference",
        "aligned": "100.00",
        "probability_component_0": "",
        "probability_component_1": "",
        "probability_component_2": "",
        "probability_main": "",
        "removed": "false",
    }
    sample_c = next(row for row in rows if row["sequence"] == "sample_c")
    assert sample_c["aligned"] == "20.00"
    assert sample_c["removed"] == "true"
    assert sample_c["probability_main"] != ""


def test_filter_alignment_by_aligned_percentage_keeps_all_samples_when_no_outliers(tmp_path):
    aln = tmp_path / "core.full.aln"
    aln.write_text(
        ">reference\nAAAAAA\n"
        ">sample_a\nAAAAAA\n"
        ">sample_b\nAAA-AAA\n"
        ">sample_c\nAA----\n"
    )
    aligned_tsv = tmp_path / "core.aligned.tsv"
    aligned_tsv.write_text(
        "sequence\taligned\n"
        "reference\t100.00\n"
        "sample_a\t100.00\n"
        "sample_b\t99.00\n"
        "sample_c\t99.50\n"
    )
    filtered_aln = tmp_path / "core.filtered.aln"

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.1,
    )
    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference", "sample_a", "sample_b", "sample_c"]
    with aligned_tsv.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sequence"] for row in rows] == ["reference", "sample_a", "sample_b", "sample_c"]
    assert all(row["removed"] == "false" for row in rows)
    assert all(row["probability_main"] == "" for row in rows)


def test_filter_alignment_by_aligned_percentage_keeps_all_samples_when_too_few_samples(tmp_path):
    aln = tmp_path / "core.full.aln"
    aln.write_text(
        ">reference\nAAAAAA\n"
        ">sample_a\nAAAAAA\n"
        ">sample_b\nAAA-AAA\n"
    )
    aligned_tsv = tmp_path / "core.aligned.tsv"
    aligned_tsv.write_text(
        "sequence\taligned\n"
        "reference\t100.00\n"
        "sample_a\t100.00\n"
        "sample_b\t85.00\n"
    )
    filtered_aln = tmp_path / "core.filtered.aln"

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference", "sample_a", "sample_b"]
    with aligned_tsv.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sequence"] for row in rows] == ["reference", "sample_a", "sample_b"]
    assert all(row["removed"] == "false" for row in rows)
    assert all(row["probability_main"] == "" for row in rows)

def test_filter_alignment_by_aligned_percentage_filters_high_and_low_outliers(tmp_path):
    aln = tmp_path / "core.full.aln"
    aln.write_text(
        ">reference\nAAAAAA\n"
        ">sample_a\nAAAAAA\n"
        ">sample_b\nAAA-AAA\n"
        ">sample_c\nAA----\n"
        ">sample_d\n------\n"
        ">sample_e\nAAAAAA\n"
        ">sample_f\nAAAAAA\n"
    )
    aligned_tsv = tmp_path / "core.aligned.tsv"
    aligned_tsv.write_text(
        "sequence\taligned\n"
        "reference\t100.00\n"
        "sample_a\t99.00\n"
        "sample_b\t85.00\n"
        "sample_c\t83.00\n"
        "sample_d\t80.00\n"
        "sample_e\t10.00\n"
        "sample_f\t35.00\n"
    )
    filtered_aln = tmp_path / "core.filtered.aln"

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference", "sample_b", "sample_c", "sample_d"]
    with aligned_tsv.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    removed_rows = {row["sequence"]: row for row in rows if row["removed"] == "true"}
    assert set(removed_rows) == {"sample_a", "sample_e", "sample_f"}
    assert removed_rows["sample_a"]["probability_main"] != ""
    assert removed_rows["sample_a"]["probability_component_2"] != ""
    kept_rows = {row["sequence"]: row for row in rows if row["removed"] == "false"}
    assert float(kept_rows["sample_b"]["probability_main"]) >= 0.50
    assert float(removed_rows["sample_a"]["probability_main"]) < 0.50

def test_filter_alignment_by_aligned_percentage_keeps_samples_when_one_tail_outlier(tmp_path):
    aln = tmp_path / "core.full.aln"
    aln.write_text(
        ">reference\nAAAAAA\n"
        ">sample_a\nAAAAAA\n"
        ">sample_b\nAAA-AAA\n"
        ">sample_c\nAA----\n"
        ">sample_d\n------\n"
        ">sample_e\nAAAAAA\n"
    )
    aligned_tsv = tmp_path / "core.aligned.tsv"
    aligned_tsv.write_text(
        "sequence\taligned\n"
        "reference\t100.00\n"
        "sample_a\t99.00\n"
        "sample_b\t99.00\n"
        "sample_c\t99.00\n"
        "sample_d\t60.00\n"
        "sample_e\t30.00\n"
    )
    filtered_aln = tmp_path / "core.filtered.aln"

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference", "sample_a", "sample_b", "sample_c"]

def test_filter_alignment_by_aligned_percentage_warns_when_samples_removed(tmp_path, monkeypatch):
    aln = tmp_path / "core.full.aln"
    aln.write_text(
        ">reference\nAAAAAA\n"
        ">sample_a\nAAAAAA\n"
        ">sample_b\nAAA-AAA\n"
        ">sample_c\nAA----\n"
    )
    aligned_tsv = tmp_path / "core.aligned.tsv"
    aligned_tsv.write_text(
        "sequence\taligned\n"
        "reference\t100.00\n"
        "sample_a\t100.00\n"
        "sample_b\t85.00\n"
        "sample_c\t20.00\n"
    )
    filtered_aln = tmp_path / "core.filtered.aln"
    messages: list[str] = []

    monkeypatch.setattr(logger, "warning", messages.append)

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    assert messages == [
        f"Filtered out 1 sample(s) from alignment {aln} using inclusion threshold 0.50: sample_c"
    ]


def test_core_pipeline_soft_core_uses_filtered_alignment(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "reference.fa").write_text(">ref\nAAAA\n")
    (reference / "reference.fa.fai").write_text("ref\t4\t0\t4\t5\n")
    (reference / "reference.dict").write_text("ref\t4\n")
    (reference / "reference.gff").write_text("##gff-version 3\n")
    (reference / "metadata.json").write_text('{"prefix": "reference"}')

    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()

    pipeline = CorePipelineBuilder(
        snippy_dirs=[snippy_dir],
        reference=reference,
        prefix="core",
    ).build()

    soft_core_stage = next(stage for stage in pipeline.stages if isinstance(stage, SoftCoreFilter))
    filter_stage = next(stage for stage in pipeline.stages if isinstance(stage, FilterAlignmentByAlignedPercentage))

    assert soft_core_stage.aln == filter_stage.output.filtered_aln


def test_distle_distance_matrix_uses_phylip_output(tmp_path):
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
        "phylip",
        str(tmp_path / "core.full.aln"),
        str(tmp_path / "core.full.phylip"),
    ]
    assert stage.output.phylip == tmp_path / "core.full.phylip"


def test_core_pipeline_adds_distle_for_full_and_soft_core_alignments(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "reference.fa").write_text(">ref\nAAAA\n")
    (reference / "reference.fa.fai").write_text("ref\t4\t0\t4\t5\n")
    (reference / "reference.dict").write_text("ref\t4\n")
    (reference / "reference.gff").write_text("##gff-version 3\n")
    (reference / "metadata.json").write_text('{"prefix": "reference"}')

    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()

    pipeline = CorePipelineBuilder(
        snippy_dirs=[snippy_dir],
        reference=reference,
        prefix="core",
    ).build()

    distle_stages = [stage for stage in pipeline.stages if isinstance(stage, DistleDistanceMatrix)]
    combine_stage = next(stage for stage in pipeline.stages if isinstance(stage, CombineFastaFile))
    soft_core_stage = next(stage for stage in pipeline.stages if isinstance(stage, SoftCoreFilter))

    assert len(distle_stages) == 2
    assert distle_stages[0].aln == combine_stage.output.aln
    assert distle_stages[0].output.phylip == Path("core.full.phylip")
    assert distle_stages[1].aln == soft_core_stage.output.soft_core
    assert distle_stages[1].output.phylip == Path("core.095.phylip")
