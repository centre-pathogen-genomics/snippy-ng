from io import StringIO

import pytest

from snippy_ng.utils.cnv import (
    CNVError,
    CoverageRow,
    FeatureRow,
    copy_number_variation,
    estimate_copy_numbers,
    estimate_feature_copy_numbers,
    parse_known_single_copy_region,
    parse_samtools_coverage,
    parse_samtools_depth,
    parse_gff_features,
    SingleCopyRegion,
    write_feature_cnv_table,
    write_cnv_table,
)


def test_parse_samtools_coverage_reads_mean_depth_rows():
    rows = parse_samtools_coverage(
        [
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq",
            "chr1\t1\t1000\t100\t1000\t100\t30.5\t40\t60",
            "plasmid\t1\t100\t40\t100\t100\t92\t40\t60",
        ]
    )

    assert rows == [
        CoverageRow(contig_id="chr1", length=1000, read_depth=30.5),
        CoverageRow(contig_id="plasmid", length=100, read_depth=92.0),
    ]


def test_estimate_copy_numbers_uses_largest_contig_as_copy_one_baseline():
    cnv_rows = estimate_copy_numbers(
        [
            CoverageRow(contig_id="chr1", length=1000, read_depth=30.0),
            CoverageRow(contig_id="chr2_dup", length=800, read_depth=61.0),
            CoverageRow(contig_id="plasmid", length=100, read_depth=136.0),
            CoverageRow(contig_id="missing", length=50, read_depth=0.0),
        ]
    )

    assert [(row.contig_id, row.copy_number) for row in cnv_rows] == [
        ("chr1", 1),
        ("chr2_dup", 2),
        ("plasmid", 5),
        ("missing", 0),
    ]


def test_estimate_copy_numbers_rejects_zero_depth_baseline():
    with pytest.raises(CNVError, match="Largest contig 'chr1' has zero depth"):
        estimate_copy_numbers(
            [
                CoverageRow(contig_id="chr1", length=1000, read_depth=0.0),
                CoverageRow(contig_id="plasmid", length=100, read_depth=20.0),
            ]
        )


def test_write_cnv_table_outputs_requested_columns():
    rows = estimate_copy_numbers(
        [
            CoverageRow(contig_id="chr1", length=1000, read_depth=30.0),
            CoverageRow(contig_id="plasmid", length=100, read_depth=90.0),
        ]
    )
    output = StringIO()

    write_cnv_table(rows, output)

    assert output.getvalue() == (
        "contig_id\tread_depth\tcopy_number\n"
        "chr1\t30\t1\n"
        "plasmid\t90\t3\n"
    )


def test_estimate_copy_numbers_accepts_known_baseline_depth():
    cnv_rows = estimate_copy_numbers(
        [
            CoverageRow(contig_id="chr1", length=1000, read_depth=60.0),
            CoverageRow(contig_id="plasmid", length=100, read_depth=120.0),
        ],
        baseline=30.0,
    )

    assert [(row.contig_id, row.copy_number) for row in cnv_rows] == [
        ("chr1", 2),
        ("plasmid", 4),
    ]


def test_parse_known_single_copy_region_accepts_largest_contig_coordinate_form():
    assert parse_known_single_copy_region("4518,5000") == SingleCopyRegion(
        start=4518,
        end=5000,
    )


def test_parse_known_single_copy_region_accepts_named_region_form():
    assert parse_known_single_copy_region("chr1:4518-5000") == SingleCopyRegion(
        contig_id="chr1",
        start=4518,
        end=5000,
    )


def test_parse_gff_features_extracts_selected_feature_type(tmp_path):
    gff = tmp_path / "ref.gff"
    gff.write_text(
        "chr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1\n"
        "chr1\t.\tCDS\t2\t9\t.\t+\t0\tID=cds1;gene=abc\n"
        "chr2\t.\tCDS\t5\t7\t.\t-\t0\tName=cds2\n"
    )

    assert parse_gff_features(gff, feature_type="CDS") == [
        FeatureRow(feature_id="cds1", contig_id="chr1", start=2, end=9),
        FeatureRow(feature_id="cds2", contig_id="chr2", start=5, end=7),
    ]


def test_parse_samtools_depth_and_estimate_feature_copy_numbers_uses_median_depth():
    depths = parse_samtools_depth(
        [
            "chr1\t1\t30",
            "chr1\t2\t31",
            "chr1\t3\t1000",
            "chr1\t4\t29",
            "chr1\t5\t30",
            "plasmid\t1\t90",
            "plasmid\t2\t91",
            "plasmid\t3\t89",
        ]
    )
    rows = estimate_feature_copy_numbers(
        [
            FeatureRow(feature_id="cds1", contig_id="chr1", start=1, end=5),
            FeatureRow(feature_id="cds2", contig_id="plasmid", start=1, end=3),
            FeatureRow(feature_id="missing", contig_id="chr1", start=10, end=11),
        ],
        depths,
        baseline=30.0,
    )

    assert [(row.feature_id, row.read_depth, row.copy_number) for row in rows] == [
        ("cds1", 30.0, 1),
        ("cds2", 90.0, 3),
        ("missing", 0.0, 0),
    ]


def test_write_feature_cnv_table_outputs_requested_columns():
    rows = estimate_feature_copy_numbers(
        [FeatureRow(feature_id="cds1", contig_id="chr1", start=1, end=3)],
        {"chr1": {1: 30, 2: 31, 3: 32}},
        baseline=30.0,
    )
    output = StringIO()

    write_feature_cnv_table(rows, output)

    assert output.getvalue() == (
        "feature_id\tcontig_id\tstart\tend\tread_depth\tcopy_number\n"
        "cds1\tchr1\t1\t3\t31\t1\n"
    )


def test_copy_number_variation_returns_contig_rows(monkeypatch, tmp_path):
    alignment = tmp_path / "sample.cram"
    alignment.write_text("cram")

    monkeypatch.setattr(
        "snippy_ng.utils.cnv.samtools_coverage",
        lambda alignment: (
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            "chr1\t1\t1000\t100\t1000\t100\t30\t40\t60\n"
            "plasmid\t1\t100\t100\t100\t100\t90\t40\t60\n"
        ),
    )

    rows = copy_number_variation(alignment)

    assert [(row.contig_id, row.read_depth, row.copy_number) for row in rows] == [
        ("chr1", 30.0, 1),
        ("plasmid", 90.0, 3),
    ]


def test_copy_number_variation_returns_feature_rows(monkeypatch, tmp_path):
    alignment = tmp_path / "sample.cram"
    alignment.write_text("cram")
    gff = tmp_path / "reference.gff"
    gff.write_text("plasmid\t.\tCDS\t1\t3\t.\t+\t0\tID=cds1\n")

    monkeypatch.setattr(
        "snippy_ng.utils.cnv.samtools_coverage",
        lambda alignment: (
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            "chr1\t1\t1000\t100\t1000\t100\t30\t40\t60\n"
            "plasmid\t1\t100\t100\t100\t100\t90\t40\t60\n"
        ),
    )
    monkeypatch.setattr(
        "snippy_ng.utils.cnv.samtools_depth",
        lambda alignment, bed: (
            "plasmid\t1\t89\n"
            "plasmid\t2\t90\n"
            "plasmid\t3\t91\n"
        ),
    )

    rows = copy_number_variation(alignment, gff=gff, feature_type="CDS")

    assert [(row.feature_id, row.read_depth, row.copy_number) for row in rows] == [
        ("cds1", 90.0, 3),
    ]
