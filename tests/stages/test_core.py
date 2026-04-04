import csv

from Bio import SeqIO

from snippy_ng.pipelines.core import CorePipelineBuilder
from snippy_ng.stages.core import FilterAlignmentByAlignedPercentage, SoftCoreFilter


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
        "probability_main": "",
        "removed": "false",
    }
    sample_c = next(row for row in rows if row["sequence"] == "sample_c")
    assert sample_c["aligned"] == "20.00"
    assert sample_c["removed"] == "true"
    assert sample_c["probability_main"] != ""


def test_filter_alignment_by_aligned_percentage_keeps_all_when_too_few_samples(tmp_path):
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
