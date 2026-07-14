import csv

import numpy as np
import pytest
from Bio import SeqIO
from pydantic import ValidationError

import snippy_ng.stages.alignment_filter as filter_module
from snippy_ng.exceptions import MSAValidationError
from snippy_ng.logging import logger
from snippy_ng.stages.alignment_filter import FilterAlignmentByAlignedPercentage


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


def _filter_stats_path(filtered_aln):
    return filtered_aln.with_suffix(".stats.tsv")


def test_filter_alignment_retains_components_above_largest_gap(tmp_path, monkeypatch):
    aln, aligned_tsv, filtered_aln = _write_filter_files(
        tmp_path,
        [20, 21, 93, 94, 94, 95, 95, 96, 99, 100],
    )
    original_stats = aligned_tsv.read_text()

    def fake_select_gmm(cls, values, **kwargs):
        responsibilities = np.array(
            [[1.0, 0.0, 0.0]] * 2
            + [[0.0, 1.0, 0.0]] * 6
            + [[0.0, 0.0, 1.0]] * 2
        )
        return (
            3,
            np.array([0.2, 0.6, 0.2]),
            np.array([20.5, 94.5, 99.5]),
            np.array([0.25, 1.0, 0.25]),
            responsibilities,
            {1: 100.0, 2: 70.0, 3: 50.0},
        )

    monkeypatch.setattr(
        FilterAlignmentByAlignedPercentage,
        "_select_gmm",
        classmethod(fake_select_gmm),
    )

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln=aln,
        alignment_stats=aligned_tsv,
        filtered_aln=filtered_aln,
        inclusion_threshold=0.50,
    )

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference"] + [f"sample_{index}" for index in range(2, 10)]
    assert aligned_tsv.read_text() == original_stats

    with _filter_stats_path(filtered_aln).open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["filter_reason"] == "reference"
    assert rows[0]["selected_components"] == "3"
    assert rows[0]["retained_components"] == "1,2"
    assert rows[1]["filter_reason"] == "removed_by_gmm"
    assert rows[-1]["filter_reason"] == "retained_by_gmm"
    assert rows[-1]["probability_component_2"] == "1.0000"
    assert rows[-1]["probability_retained"] == "1.0000"
    assert rows[-1]["probability_main"] == rows[-1]["probability_retained"]


def test_filter_alignment_keeps_both_components_when_separation_is_small(tmp_path, monkeypatch):
    aln, aligned_tsv, filtered_aln = _write_filter_files(
        tmp_path,
        [89, 94, 94, 95, 95, 96, 96, 97, 99, 100],
    )

    def fake_select_gmm(cls, values, **kwargs):
        return (
            2,
            np.array([0.8, 0.2]),
            np.array([94.0, 100.0]),
            np.array([2.0, 0.5]),
            np.array([[1.0, 0.0]] * 8 + [[0.0, 1.0]] * 2),
            {1: 80.0, 2: 60.0},
        )

    monkeypatch.setattr(
        FilterAlignmentByAlignedPercentage,
        "_select_gmm",
        classmethod(fake_select_gmm),
    )

    FilterAlignmentByAlignedPercentage.filter_alignment(
        aln,
        aligned_tsv,
        filtered_aln,
    )

    kept_ids = [record.id for record in SeqIO.parse(str(filtered_aln), "fasta")]
    assert kept_ids == ["reference"] + [f"sample_{index}" for index in range(10)]
    with _filter_stats_path(filtered_aln).open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["filter_reason"] == "low_component_separation" for row in rows[1:])
    assert rows[0]["retained_components"] == "0,1"
    assert all(row["probability_retained"] == "1.0000" for row in rows[1:])


@pytest.mark.parametrize(
    ("bics", "expected_components"),
    [
        ({1: 100.0, 2: 93.0, 3: 88.0}, 2),
        ({1: 100.0, 2: 93.0, 3: 86.0}, 3),
    ],
)
def test_select_gmm_requires_bic_improvement_for_complexity(monkeypatch, bics, expected_components):
    def fake_fit_gmm(cls, values, n_components=2):
        return (
            np.full(n_components, 1.0 / n_components),
            np.arange(n_components, dtype=float),
            np.ones(n_components),
            np.full((len(values), n_components), 1.0 / n_components),
            bics[n_components],
        )

    monkeypatch.setattr(
        FilterAlignmentByAlignedPercentage,
        "_fit_gmm",
        classmethod(fake_fit_gmm),
    )

    selected, *_ = FilterAlignmentByAlignedPercentage._select_gmm(
        list(range(15)),
        bic_improvement=6.0,
    )

    assert selected == expected_components


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
    with _filter_stats_path(filtered_aln).open("r", newline="") as handle:
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
    with _filter_stats_path(filtered_aln).open("r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sequence"] for row in rows] == ["reference", "sample_a", "sample_b"]
    assert all(row["removed"] == "false" for row in rows)
    assert all(row["probability_main"] == "" for row in rows)
    assert all(row["filter_reason"] == "too_few_samples" for row in rows[1:])


def test_filter_alignment_by_aligned_percentage_rejects_missing_sample_stats(tmp_path):
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
    with pytest.raises(ValueError, match="Missing aligned percentage for sequence 'sample_b'"):
        FilterAlignmentByAlignedPercentage.filter_alignment(
            aln=aln,
            alignment_stats=aligned_tsv,
            filtered_aln=filtered_aln,
            inclusion_threshold=0.50,
        )


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
    with _filter_stats_path(filtered_aln).open("r", newline="") as handle:
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
        ("minimum_cluster_separation", -0.01),
        ("minimum_cluster_separation", 100.01),
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
    class NonFiniteGaussianMixture:
        def __init__(self, **kwargs):
            pass

        def fit(self, values):
            self.converged_ = True
            self.weights_ = np.array([0.5, 0.5])
            self.means_ = np.array([[0.0], [100.0]])
            self.covariances_ = np.array([[[1.0]], [[1.0]]])
            return self

        def predict_proba(self, values):
            return np.full((len(values), 2), np.nan)

        def bic(self, values):
            return 0.0

    monkeypatch.setattr(filter_module, "GaussianMixture", NonFiniteGaussianMixture)

    with pytest.raises(ValueError, match="non-finite responsibilities"):
        FilterAlignmentByAlignedPercentage._fit_gmm([0, 1, 99, 100])


def test_fit_gmm_rejects_near_empty_component():
    with pytest.raises(ValueError, match="near-empty component"):
        FilterAlignmentByAlignedPercentage._fit_gmm([0] * 99 + [100])
