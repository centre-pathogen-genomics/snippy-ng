from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.masks import DepthBedsFromBam, ApplyDepthMaskToFasta, MaskMixedSites, MinDepthThresholdFromBam
from snippy_ng.stages.vcf import AddDeletionsToVCF


def _write_minimal_vcf(vcf_path: Path) -> None:
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Wildtype,length=500000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n"
    )


def _write_reference_fasta(reference_path: Path) -> None:
    reference_path.write_text(
        ">Wildtype\n"
        + ("A" * 600000)
        + "\n"
    )


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

    AddDeletionsToVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1

    chrom, pos, _id, ref, alt, _qual, _filt, info, fmt, sample = rows[0]
    assert chrom == "Wildtype"
    assert pos == "103379"
    assert alt == "<DEL>"
    assert "ZERODEPTH" in info
    assert "END=103380" in info
    assert "SVLEN=-1" in info
    assert fmt == "GT"
    assert sample == "1/1"


def test_add_deletions_to_vcf_contig_start_interval_uses_anchor_semantics(tmp_path: Path):
    input_vcf = tmp_path / "input.vcf"
    bed = tmp_path / "zero.bed"
    out_vcf = tmp_path / "out.vcf"
    reference = tmp_path / "reference.fa"

    _write_minimal_vcf(input_vcf)
    _write_reference_fasta(reference)
    bed.write_text("Wildtype\t0\t20\n")

    AddDeletionsToVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

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

    AddDeletionsToVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1
    chrom, pos, _id, ref, alt, *_ = rows[0]
    assert chrom == "Wildtype"
    assert pos == "103379"
    assert ref == "GA"
    assert alt == "G"


def test_add_deletions_to_vcf_skips_zero_depth_del_if_overlapping_existing_variant(tmp_path: Path):
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

    AddDeletionsToVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1
    chrom, pos, _id, ref, alt, *_ = rows[0]
    assert chrom == "Wildtype"
    assert pos == "565"
    assert ref == "A"
    assert alt == "G"


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
    MinDepthThresholdFromBam.write_min_depth_threshold(depth_txt, out_disabled, min_depth=0)
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
        "awk", '$3>0 && $3<47 {print $1"\\t"($2-1)"\\t"$2}'
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
        "samtools", "depth", "-aa", "sample.bam",
    ]
    assert commands[0].processes[1].command == [
        "awk", '$3==0 {print $1"\\t"($2-1)"\\t"$2}'
    ]
    assert commands[0].processes[2].command == ["bedtools", "merge", "-i", "-"]
    assert commands[0].output_file == Path("sample.zerodepth.bed")

    assert len(commands[1].processes) == 3
    assert commands[1].processes[0].command == [
        "samtools", "depth", "-aa", "sample.bam",
    ]
    assert commands[1].processes[1].command == [
        "awk", '$3>0 && $3<10 {print $1"\\t"($2-1)"\\t"$2}'
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
        "samtools", "depth", "-aa", "sample.bam",
    ], "Default quality thresholds should not be included in Zero-depth command"
    assert commands[1].processes[0].command == [
        "samtools", "depth", "-aa", "sample.bam",
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
        "bedtools", "maskfasta",
        "-fi", "sample.consensus.fasta",
        "-bed", "sample.mindepth.bed",
        "-fo", "sample.mindepth_masked.fasta",
        "-fullHeader",
        "-mc", "N",
    ]


def test_mixed_site_mask_uses_maskfasta_with_lowercase_n():
    stage = MaskMixedSites(
        fasta=Path("sample.consensus.fasta"),
        vcf=Path("sample.all.vcf"),
        prefix="sample",
    )

    commands = stage.create_commands(Context())

    assert commands[1].command == [
        "bedtools", "maskfasta",
        "-fi", "sample.consensus.fasta",
        "-bed", "sample.mixed_sites.bed",
        "-fo", "sample.mixed_masked.fasta",
        "-fullHeader",
        "-mc", "n",
    ]
