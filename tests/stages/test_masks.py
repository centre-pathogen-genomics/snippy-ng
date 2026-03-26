from pathlib import Path

from snippy_ng.stages.vcf import AddDeletionstoVCF


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

    AddDeletionstoVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

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

    AddDeletionstoVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

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

    AddDeletionstoVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

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

    AddDeletionstoVCF._merge_zero_depth_deletions_into_vcf(input_vcf, bed, out_vcf, reference)

    rows = _get_variant_rows(out_vcf)
    assert len(rows) == 1
    chrom, pos, _id, ref, alt, *_ = rows[0]
    assert chrom == "Wildtype"
    assert pos == "565"
    assert ref == "A"
    assert alt == "G"
