import csv
from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO
from pydantic import ValidationError

from snippy_ng.logging import logger
from snippy_ng.pipelines.core import CorePipelineBuilder
from snippy_ng.stages.core import CombineFastaFile, DistleDistanceMatrix, FilterAlignmentByAlignedPercentage, MSAValidationError, SoftCoreFilter
from snippy_ng.context import Context


def _write_filter_files(tmp_path, sample_values):
    aln = tmp_path / "core.full.aln"
    aln.write_text(">reference\nAAAAAA\n" + "".join(
        f">sample_{index}\nAAAAAA\n" for index in range(len(sample_values))
    ))
    aligned_tsv = tmp_path / "core.aligned.tsv"
    aligned_tsv.write_text("sequence\taligned\nreference\t100.00\n" + "".join(
        f"sample_{index}\t{value:.2f}\n" for index, value in enumerate(sample_values)
    ))
    filtered_aln = tmp_path / "core.filtered.aln"
    return aln, aligned_tsv, filtered_aln


def test_filter_alignment_uses_highest_mean_component(tmp_path, monkeypatch):
    aln, aligned_tsv, filtered_aln = _write_filter_files(
        tmp_path,
        [20, 21, 22, 23, 24, 25, 26, 27, 90, 91],
    )

    def fake_fit_gmm(cls, values, max_iter=100, tol=1e-6, n_components=2):
        assert n_components == 2
        responsibilities = np.array([[1.0, 0.0]] * 8 + [[0.0, 1.0]] * 2)
        return (
            np.array([0.8, 0.2]),
            np.array([23.5, 90.5]),
            np.array([5.25, 0.25]),
            responsibilities,
        )

    monkeypatch.setattr(
        FilterAlignmentByAlignedPercentage,
        "_fit_gmm",
        classmethod(fake_fit_gmm),
    )

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference", "sample_8", "sample_9"]

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
        "filter_reason": "reference",
    }
    assert rows[1]["filter_reason"] == "removed_by_gmm"
    assert rows[-1]["filter_reason"] == "retained_by_gmm"
    assert rows[-1]["probability_component_2"] == ""


def test_filter_alignment_by_aligned_percentage_keeps_all_samples_when_no_outliers(tmp_path):
    values = [99.0, 99.1, 99.2, 99.3, 99.4, 99.5, 99.6, 99.7, 99.8, 100.0]
    aln, aligned_tsv, filtered_aln = _write_filter_files(tmp_path, values)

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.1,
    )
    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference"] + [f"sample_{index}" for index in range(10)]
    with aligned_tsv.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["removed"] == "false" for row in rows)
    assert all(row["probability_main"] == "" for row in rows)
    assert rows[0]["filter_reason"] == "reference"
    assert all(row["filter_reason"] == "low_spread" for row in rows[1:])


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
    assert all(row["filter_reason"] == "too_few_samples" for row in rows[1:])


def test_filter_alignment_by_aligned_percentage_treats_missing_sample_stats_as_zero(tmp_path, monkeypatch):
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

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference", "sample_a", "sample_b"]
    with aligned_tsv.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sample_b = next(row for row in rows if row["sequence"] == "sample_b")
    assert sample_b["aligned"] == "0.00"
    assert messages == [
        f"Missing aligned percentage for sequence 'sample_b' in {aligned_tsv}; treating it as 0.00"
    ]


def test_filter_alignment_by_aligned_percentage_does_not_fit_small_dataset(tmp_path):
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
    assert kept_ids == ["reference", "sample_a", "sample_b", "sample_c", "sample_d", "sample_e", "sample_f"]
    with aligned_tsv.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["removed"] == "false" for row in rows)
    assert all(row["filter_reason"] == "too_few_samples" for row in rows[1:])

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
    assert kept_ids == ["reference", "sample_a", "sample_b", "sample_c", "sample_d", "sample_e"]

def test_filter_alignment_by_aligned_percentage_warns_when_samples_removed(tmp_path, monkeypatch):
    aln, aligned_tsv, filtered_aln = _write_filter_files(
        tmp_path,
        [99, 98, 97, 96, 95, 94, 93, 92, 20, 10],
    )
    messages: list[str] = []

    monkeypatch.setattr(logger, "warning", messages.append)

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    assert messages == [
        f"Filtered out 2 sample(s) from alignment {aln} using inclusion threshold 0.50: sample_8, sample_9"
    ]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("inclusion_threshold", -0.01),
        ("inclusion_threshold", 1.01),
        ("identical_alignment_spread", -0.01),
        ("identical_alignment_spread", 100.01),
    ],
)
def test_filter_alignment_validates_configuration(tmp_path, field, value):
    with pytest.raises(ValidationError):
        FilterAlignmentByAlignedPercentage(
            aln=tmp_path / "core.full.aln",
            alignment_stats=tmp_path / "core.aligned.tsv",
            prefix="core",
            **{field: value},
        )


@pytest.mark.parametrize("aligned", ["nan", "-0.01", "100.01"])
def test_filter_alignment_rejects_invalid_aligned_percentage(tmp_path, aligned):
    aln, aligned_tsv, filtered_aln = _write_filter_files(tmp_path, [99])
    aligned_tsv.write_text(
        f"sequence\taligned\nreference\t100.00\nsample_0\t{aligned}\n"
    )

    with pytest.raises(ValueError, match="Invalid aligned percentage for 'sample_0'"):
        FilterAlignmentByAlignedPercentage.filter_alignment(
            aln, aligned_tsv, filtered_aln
        )


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("", "Alignment statistics file has no header"),
        ("sequence\nreference\n", "Alignment statistics file is missing columns"),
    ],
)
def test_filter_alignment_validates_statistics_header(tmp_path, contents, message):
    aln, aligned_tsv, filtered_aln = _write_filter_files(tmp_path, [99])
    aligned_tsv.write_text(contents)

    with pytest.raises(ValueError, match=message):
        FilterAlignmentByAlignedPercentage.filter_alignment(
            aln, aligned_tsv, filtered_aln
        )


def test_filter_alignment_rejects_duplicate_statistics_sequences(tmp_path):
    aln, aligned_tsv, filtered_aln = _write_filter_files(tmp_path, [99])
    aligned_tsv.write_text(
        "sequence\taligned\nreference\t100.00\nsample_0\t99.00\nsample_0\t98.00\n"
    )

    with pytest.raises(ValueError, match="Duplicate sequence 'sample_0'"):
        FilterAlignmentByAlignedPercentage.filter_alignment(
            aln, aligned_tsv, filtered_aln
        )


def test_filter_alignment_rejects_duplicate_alignment_ids(tmp_path):
    aln, aligned_tsv, filtered_aln = _write_filter_files(tmp_path, [99])
    aln.write_text(">reference\nAAAAAA\n>sample_0\nAAAAAA\n>sample_0\nAAAAAA\n")

    with pytest.raises(MSAValidationError, match="Duplicate sequence IDs"):
        FilterAlignmentByAlignedPercentage.filter_alignment(
            aln, aligned_tsv, filtered_aln
        )


def test_fit_gmm_rejects_non_finite_responsibilities(monkeypatch):
    monkeypatch.setattr(np, "exp", lambda values: np.full_like(values, np.nan))

    with pytest.raises(ValueError, match="non-finite responsibilities"):
        FilterAlignmentByAlignedPercentage._fit_gmm([0, 1, 99, 100], max_iter=1)


def test_fit_gmm_rejects_near_empty_component():
    with pytest.raises(ValueError, match="near-empty component"):
        FilterAlignmentByAlignedPercentage._fit_gmm([0] * 99 + [100])


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
    assert distle_stages[0].output.phylip == Path("core.full.distance.tsv")
    assert distle_stages[1].aln == soft_core_stage.output.soft_core
    assert distle_stages[1].output.phylip == Path("core.095.distance.tsv")
