import re

from click.testing import CliRunner

from snippy_ng.cli.samcrop_cli import samcrop
from snippy_ng.utils.samcrop import samcrop_filter_lines, parse_bed_lines


HEADER = [
    "@HD\tVN:1.6\tSO:coordinate\n",
    "@SQ\tSN:ref\tLN:1000\n",
]


def make_sam(qname="r1", pos=101, cigar="100M", seq=None, tags="MD:Z:100\tNM:i:0"):
    if seq is None:
        query_len = 0
        for length, op in __import__("re").findall(r"(\d+)([MIS=X])", cigar):
            query_len += int(length)
        seq = "A" * query_len
    fields = [
        qname,
        "0",
        "ref",
        str(pos),
        "60",
        cigar,
        "*",
        "0",
        "0",
        seq,
        "I" * len(seq),
    ]
    if tags:
        fields.extend(tags.split("\t"))
    return "\t".join(fields) + "\n"


def body(lines):
    return [line for line in lines if not line.startswith("@")]


def test_samcrop_crops_100m_left_right_and_both_sides():
    intervals = parse_bed_lines(["ref\t120\t180\n", "ref\t120\t220\n", "ref\t80\t180\n"])

    both = body(samcrop_filter_lines(HEADER + [make_sam("both")], {"ref": [(120, 180)]}))[0].split("\t")
    left = body(samcrop_filter_lines(HEADER + [make_sam("left")], {"ref": [(120, 220)]}))[0].split("\t")
    right = body(samcrop_filter_lines(HEADER + [make_sam("right")], {"ref": [(80, 180)]}))[0].split("\t")

    assert intervals["ref"] == [(80, 220)]
    assert both[3] == "121"
    assert both[5] == "20H60M20H"
    assert len(both[9]) == 60
    assert left[3] == "121"
    assert left[5] == "20H80M"
    assert right[3] == "101"
    assert right[5] == "80M20H"


def test_samcrop_splits_read_across_disjoint_windows():
    out = body(samcrop_filter_lines(
        HEADER + [make_sam("split")],
        {"ref": [(120, 140), (170, 190)]},
    ))
    first = out[0].split("\t")
    second = out[1].split("\t")

    assert len(out) == 2
    assert first[0] == "split"
    assert first[3] == "121"
    assert first[5] == "20H20M60H"
    assert len(first[9]) == 20
    assert second[0] == "split"
    assert second[3] == "171"
    assert second[5] == "70H20M10H"
    assert len(second[9]) == 20


def test_samcrop_combines_existing_soft_and_hard_clips():
    record = make_sam(pos=101, cigar="5H5S90M")
    out = body(samcrop_filter_lines(HEADER + [record], {"ref": [(120, 170)]}))[0].split("\t")

    assert out[3] == "121"
    assert out[5] == "30H50M20H"
    assert len(out[9]) == 50


def test_samcrop_keeps_insertions_at_left_boundary():
    record = make_sam(pos=101, cigar="50M5I50M")
    out = body(samcrop_filter_lines(HEADER + [record], {"ref": [(150, 160)]}))[0].split("\t")

    assert out[3] == "151"
    assert out[5] == "50H5I10M40H"
    assert len(out[9]) == 15


def test_samcrop_trims_deletions_crossing_boundaries_without_zero_ops():
    record = make_sam(pos=101, cigar="50M10D50M")
    out = body(samcrop_filter_lines(HEADER + [record], {"ref": [(145, 165)]}))[0].split("\t")

    assert out[3] == "146"
    assert out[5] == "45H5M10D5M45H"
    assert all(int(length) > 0 for length, _op in re.findall(r"(\d+)([A-Z=])", out[5]))


def test_samcrop_drops_reads_with_no_overlap_or_no_retained_query():
    no_overlap = make_sam("no_overlap", pos=101, cigar="100M")
    deletion_only = make_sam("deletion_only", pos=101, cigar="50M10D50M")

    assert body(samcrop_filter_lines(HEADER + [no_overlap], {"ref": [(300, 400)]})) == []
    assert body(samcrop_filter_lines(HEADER + [deletion_only], {"ref": [(152, 155)]})) == []


def test_samcrop_removes_stale_md_and_nm_tags():
    record = make_sam(pos=101, cigar="100M", tags="MD:Z:100\tNM:i:0\tAS:i:42")
    out = body(samcrop_filter_lines(HEADER + [record], {"ref": [(120, 180)]}))[0].split("\t")

    assert "MD:Z:100" not in out
    assert "NM:i:0" not in out
    assert any(field.startswith("AS:i:42") for field in out)


def test_samcrop_cli_accepts_file_and_bed(tmp_path):
    sam = tmp_path / "input.sam"
    bed = tmp_path / "windows.bed"
    sam.write_text("".join(HEADER + [make_sam()]))
    bed.write_text("ref\t120\t180\n")

    result = CliRunner().invoke(samcrop, ["--bed", str(bed), str(sam)])

    assert result.exit_code == 0
    assert "20H60M20H" in result.output


def test_samcrop_cli_accepts_stdin(tmp_path):
    bed = tmp_path / "windows.bed"
    bed.write_text("ref\t120\t180\n")

    result = CliRunner().invoke(samcrop, ["--bed", str(bed)], input="".join(HEADER + [make_sam()]))

    assert result.exit_code == 0
    assert "20H60M20H" in result.output
