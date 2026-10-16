from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.masks import (
    DepthBedsFromBam,
    ApplyDepthMaskToFasta,
    MaskMixedSites,
    MinDepthThresholdFromBam,
)
from snippy_ng.stages.vcf import AddDeletionsToVCF


def _write_minimal_vcf(vcf_path: Path) -> None:
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Wildtype,length=500000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n"
    )


def _write_reference_fasta(reference_path: Path) -> None:
    reference_path.write_text(">Wildtype\n" + ("A" * 600000) + "\n")


def _get_variant_rows(vcf_path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with open(vcf_path, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            rows.append(line.rstrip("\n").split("\t"))
    return rows


def test_add_deletions_to_vcf_symbolic_len_standard_interval(tmp_path: Path):
    input_vcf = tmp_path / "input.vcf"
    bed = tmp_path / "zero.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    _write_minimal_vcf(input_vcf)
    _write_reference_fasta(reference)
    bed.write_text("Wildtype\t103379\t103380\n")

    AddDeletionsToVCF._merge_depth_deletions_into_vcf(
        input_vcf, bed, None, out_vcf, reference
    )

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1

    chrom, pos, _id, ref, alt, _qual, _filt, info, fmt, sample = rows[0]
    assert chrom == "Wildtype"
    assert pos == "103380"
    assert ref == "A"
    assert alt == "-"
    assert "ZERODEPTH" in info
    assert "END=103380" in info
    assert "SVLEN=-1" in info
    assert fmt == "GT"
    assert sample == "1/1"


def test_add_deletions_to_vcf_contig_start_interval_uses_anchor_semantics(
    tmp_path: Path,
):
    input_vcf = tmp_path / "input.vcf"
    bed = tmp_path / "zero.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    _write_minimal_vcf(input_vcf)
    _write_reference_fasta(reference)
    bed.write_text("Wildtype\t0\t20\n")

    AddDeletionsToVCF._merge_depth_deletions_into_vcf(
        input_vcf, bed, None, out_vcf, reference
    )

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1

    chrom, pos, _id, ref, alt, _qual, _filt, info, fmt, sample = rows[0]
    assert chrom == "Wildtype"
    assert pos == "1"
    assert len(ref) == 21
    assert alt == "-"
    assert "ZERODEPTH" in info
    assert "END=20" in info
    assert "SVLEN=-20" in info
    assert fmt == "GT"
    assert sample == "1/1"

    # Ensure expected INFO headers are auto-added
    header_text = out_vcf.read_text()
    assert "##INFO=<ID=SVLEN" in header_text
    assert "##INFO=<ID=TYPE" in header_text
    assert "##INFO=<ID=ZERODEPTH" in header_text


def test_add_deletions_to_vcf_skips_existing_explicit_deletion(tmp_path: Path):
    input_vcf = tmp_path / "input.vcf"
    bed = tmp_path / "zero.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Wildtype,length=500000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n"
        "Wildtype\t103379\t.\tGA\tG\t319.528\t.\t"
        "AB=0;AO=11;DP=11;QA=373;QR=0;RO=0;TYPE=INDEL\t"
        "GT:GQ:DP:RO:QR:AO:QA:GL\t1:86:11:0:0:11:373:-33.8928,-3.31133,0\n"
    )
    _write_reference_fasta(reference)
    bed.write_text("Wildtype\t103379\t103380\n")

    AddDeletionsToVCF._merge_depth_deletions_into_vcf(
        input_vcf, bed, None, out_vcf, reference
    )

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1
    chrom, pos, _id, ref, alt, *_ = rows[0]
    assert chrom == "Wildtype"
    assert pos == "103379"
    assert ref == "GA"
    assert alt == "G"


def test_add_deletions_to_vcf_skips_zero_depth_del_if_overlapping_existing_variant(
    tmp_path: Path,
):
    input_vcf = tmp_path / "input.vcf"
    bed = tmp_path / "zero.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Wildtype,length=500000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n"
        "Wildtype\t565\t.\tA\tG\t15.83\tPASS\tF\tGT:GQ:DP:AD:AF\t0/1:15:4:3,1:0.2500\n"
    )
    _write_reference_fasta(reference)
    bed.write_text("Wildtype\t560\t570\n")

    AddDeletionsToVCF._merge_depth_deletions_into_vcf(
        input_vcf, bed, None, out_vcf, reference
    )

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1
    chrom, pos, _id, ref, alt, *_ = rows[0]
    assert chrom == "Wildtype"
    assert pos == "565"
    assert ref == "A"
    assert alt == "G"


def test_add_deletions_to_vcf_adds_low_depth_blocks(tmp_path: Path):
    input_vcf = tmp_path / "input.vcf"
    zero_bed = tmp_path / "zero.bed"
    low_bed = tmp_path / "low.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    _write_minimal_vcf(input_vcf)
    _write_reference_fasta(reference)
    zero_bed.write_text("Wildtype\t100\t101\n")
    low_bed.write_text("Wildtype\t200\t300\n")

    AddDeletionsToVCF._merge_depth_deletions_into_vcf(
        input_vcf, zero_bed, low_bed, out_vcf, reference
    )

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 2

    (
        zero_chrom,
        zero_pos,
        zero_id,
        zero_ref,
        zero_alt,
        _q,
        zero_filt,
        zero_info,
        fmt,
        sample,
    ) = rows[0]
    assert zero_pos == "101"
    assert zero_ref == "A"
    assert zero_alt == "-"
    assert zero_id == "DEL_1"
    assert zero_filt == "PASS;ZERODEPTH"
    assert "ZERODEPTH" in zero_info
    assert "SVTYPE=DEL" in zero_info

    (
        low_chrom,
        low_pos,
        low_id,
        low_ref,
        low_alt,
        _q,
        low_filt,
        low_info,
        fmt,
        sample,
    ) = rows[1]
    assert low_pos == "201"
    assert low_ref == "A" * 100
    assert low_alt == "N" * 100
    assert low_id == "LOWDEPTH_1"
    assert low_filt == "LOWDEPTHSPAN"
    assert "LOWDEPTH" in low_info
    assert "END=300" in low_info
    # LOWDEPTH spans are N-mask regions, not deletions
    assert "SVTYPE=DEL" not in low_info
    assert "SVLEN" not in low_info
    assert "TYPE=INDEL" not in low_info
    assert fmt == "GT"
    assert sample == "1/1"

    # IDs are unique across both block types
    assert zero_id != low_id
    # LOWDEPTH INFO and LOWDEPTHSPAN FILTER headers are auto-added
    header_text = out_vcf.read_text()
    assert "##INFO=<ID=LOWDEPTH" in header_text
    assert "##INFO=<ID=ZERODEPTH" in header_text
    assert "##FILTER=<ID=LOWDEPTHSPAN" in header_text


def test_add_deletions_to_vcf_spans_low_depth_over_existing_variant(tmp_path: Path):
    input_vcf = tmp_path / "input.vcf"
    zero_bed = tmp_path / "zero.bed"
    low_bed = tmp_path / "low.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Wildtype,length=500000>\n"
        '##FILTER=<ID=LOWDEPTH,Description="Low depth">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n"
        "Wildtype\t565\t.\tA\tG\t15.83\tLOWDEPTH\tDP=14\tGT:GQ:DP:AD:AF\t0/1:15:4:3,1:0.2500\n"
    )
    _write_reference_fasta(reference)
    zero_bed.write_text("")
    low_bed.write_text("Wildtype\t560\t570\n")

    AddDeletionsToVCF._merge_depth_deletions_into_vcf(
        input_vcf, zero_bed, low_bed, out_vcf, reference
    )

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 2

    # the LOWDEPTH span covers the whole low-depth interval including the variant
    span_chrom, span_pos, span_id, span_ref, span_alt, *_ = rows[0]
    assert span_chrom == "Wildtype"
    assert span_pos == "561"
    assert span_ref == "A" * 10
    assert span_alt == "N" * 10
    assert span_id == "LOWDEPTH_1"

    var_chrom, var_pos, _id, var_ref, var_alt, *_ = rows[1]
    assert var_chrom == "Wildtype"
    assert var_pos == "565"
    assert var_ref == "A"
    assert var_alt == "G"


def test_add_deletions_to_vcf_low_depth_span_sorts_before_same_pos_variant(
    tmp_path: Path,
):
    """The span must precede a variant at the same POS so bcftools consensus applies the N-span."""
    input_vcf = tmp_path / "input.vcf"
    zero_bed = tmp_path / "zero.bed"
    low_bed = tmp_path / "low.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Wildtype,length=500000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n"
        "Wildtype\t561\t.\tA\tG\t15.83\tLOWDEPTH\tDP=14\tGT:GQ:DP:AD:AF\t0/1:15:4:3,1:0.2500\n"
    )
    _write_reference_fasta(reference)
    zero_bed.write_text("")
    low_bed.write_text("Wildtype\t560\t570\n")

    AddDeletionsToVCF._merge_depth_deletions_into_vcf(
        input_vcf, zero_bed, low_bed, out_vcf, reference
    )

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 2
    assert rows[0][0] == "Wildtype"
    assert rows[0][3] == "A" * 10  # span REF
    assert rows[0][4] == "N" * 10  # span ALT
    assert rows[1][1] == "561"
    assert rows[1][4] == "G"


def test_min_depth_threshold_from_bam_uses_samtools_depth_and_python():
    stage = MinDepthThresholdFromBam(
        bam=Path("sample.bam"),
        min_depth=5,
        prefix="sample",
    )

    commands = stage.create_commands(Context())

    assert len(commands) == 2
    assert commands[0].command == ["samtools", "depth", "sample.bam"]
    assert commands[0].output_file == Path("sample.mindepth.depth.txt")
    assert commands[1].func == MinDepthThresholdFromBam.write_min_depth_threshold
    assert commands[1].args == [
        Path("sample.mindepth.depth.txt"),
        Path("sample.mindepth.threshold.txt"),
        5,
    ]
    assert stage.output.threshold_txt == Path("sample.mindepth.threshold.txt")


def test_write_min_depth_threshold_median_minus_2mad(tmp_path):
    depth_txt = tmp_path / "depth.txt"
    depth_txt.write_text(
        "chr1\t1\t60\n"
        "chr1\t2\t100\n"
        "chr1\t3\t100\n"
        "chr1\t4\t100\n"
        "chr1\t5\t140\n"
        "chr1\t6\t140\n"
    )

    out = tmp_path / "threshold.txt"
    MinDepthThresholdFromBam.write_min_depth_threshold(depth_txt, out, min_depth=5)
    assert out.read_text().strip() == "60"  # median 100, MAD 20 -> 100 - 2*20

    floored = tmp_path / "floored.txt"
    MinDepthThresholdFromBam.write_min_depth_threshold(depth_txt, floored, min_depth=70)
    assert floored.read_text().strip() == "70"  # user lower bound wins


def test_write_min_depth_threshold_empty_depth_falls_back_to_min_depth(tmp_path):
    depth_txt = tmp_path / "depth.txt"
    depth_txt.write_text("")

    out = tmp_path / "threshold.txt"
    MinDepthThresholdFromBam.write_min_depth_threshold(depth_txt, out, min_depth=10)
    assert out.read_text().strip() == "10"  # nothing covered: floor applies

    out_disabled = tmp_path / "disabled.txt"
    MinDepthThresholdFromBam.write_min_depth_threshold(
        depth_txt, out_disabled, min_depth=0
    )
    assert out_disabled.read_text().strip() == "0"  # no floor: masking disabled


def test_depth_beds_from_bam_uses_dynamic_min_depth_file(tmp_path):
    threshold = tmp_path / "threshold.txt"
    threshold.write_text("47\n")

    stage = DepthBedsFromBam(
        bam=Path("sample.bam"),
        min_depth=0,
        min_depth_file=threshold,
        prefix="sample",
    )

    commands = stage.create_commands(Context())
    assert commands[1].processes[1].command == [
        "awk",
        '$3>0 && $3<47 {print $1"\\t"($2-1)"\\t"$2}',
    ]


def test_depth_mask_from_bed_uses_dynamic_min_depth_file(tmp_path):
    threshold = tmp_path / "threshold.txt"
    threshold.write_text("47\n")

    stage = ApplyDepthMaskToFasta(
        fasta=Path("sample.consensus.fasta"),
        mask_bed=Path("sample.mindepth.bed"),
        min_depth=0,
        min_depth_file=threshold,
        prefix="sample",
    )

    commands = stage.create_commands(Context())
    assert "0 < depth < 47" in commands[0].description


def test_depth_beds_from_bam_uses_samtools_depth_pipelines():
    stage = DepthBedsFromBam(
        bam=Path("sample.bam"),
        min_depth=10,
        prefix="sample",
    )

    commands = stage.create_commands(Context())

    assert len(commands) == 2
    assert len(commands[0].processes) == 3
    assert commands[0].processes[0].command == [
        "samtools",
        "depth",
        "-aa",
        "sample.bam",
    ]
    assert commands[0].processes[1].command == [
        "awk",
        '$3==0 {print $1"\\t"($2-1)"\\t"$2}',
    ]
    assert commands[0].processes[2].command == ["bedtools", "merge", "-i", "-"]
    assert commands[0].output_file == Path("sample.zerodepth.bed")

    assert len(commands[1].processes) == 3
    assert commands[1].processes[0].command == [
        "samtools",
        "depth",
        "-aa",
        "sample.bam",
    ]
    assert commands[1].processes[1].command == [
        "awk",
        '$3>0 && $3<10 {print $1"\\t"($2-1)"\\t"$2}',
    ]
    assert commands[1].processes[2].command == ["bedtools", "merge", "-i", "-"]
    assert commands[1].output_file == Path("sample.mindepth.bed")


def test_depth_beds_from_bam_can_match_variant_calling_quality_thresholds():
    stage = DepthBedsFromBam(
        bam=Path("sample.bam"),
        min_depth=10,
        prefix="sample",
    )

    commands = stage.create_commands(Context())

    assert commands[0].processes[0].command == [
        "samtools",
        "depth",
        "-aa",
        "sample.bam",
    ], "Default quality thresholds should not be included in Zero-depth command"
    assert commands[1].processes[0].command == [
        "samtools",
        "depth",
        "-aa",
        "sample.bam",
    ]


def test_depth_mask_from_bed_uses_maskfasta_with_n():
    stage = ApplyDepthMaskToFasta(
        fasta=Path("sample.consensus.fasta"),
        mask_bed=Path("sample.mindepth.bed"),
        min_depth=10,
        prefix="sample",
    )

    commands = stage.create_commands(Context())

    assert len(commands) == 1
    assert commands[0].command == [
        "bedtools",
        "maskfasta",
        "-fi",
        "sample.consensus.fasta",
        "-bed",
        "sample.mindepth.bed",
        "-fo",
        "sample.mindepth_masked.fasta",
        "-fullHeader",
        "-mc",
        "N",
    ]


def test_mixed_site_mask_uses_maskfasta_with_lowercase_n():
    stage = MaskMixedSites(
        fasta=Path("sample.consensus.fasta"),
        vcf=Path("sample.all.vcf"),
        prefix="sample",
    )

    commands = stage.create_commands(Context())

    assert commands[1].command == [
        "bedtools",
        "maskfasta",
        "-fi",
        "sample.consensus.fasta",
        "-bed",
        "sample.mixed_sites.bed",
        "-fo",
        "sample.mixed_masked.fasta",
        "-fullHeader",
        "-mc",
        "n",
    ]
