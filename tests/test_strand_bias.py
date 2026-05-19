from pathlib import Path

from snippy_ng.utils.strand_bias import (
    _format_pvalue_and_phred,
    calculate_strand_bias,
    count_strands_for_variant,
    fisher_exact_two_sided,
    annotate_vcf_strand_bias,
)


def test_fisher_exact_two_sided_returns_one_for_balanced_table():
    assert fisher_exact_two_sided(5, 5, 3, 3) == 1.0


def test_count_strands_for_snp_uses_case_and_reference_markers():
    counts = count_strands_for_variant(ref="A", alt="G", read_bases="..,,Gg")

    assert counts is not None
    assert (counts.ref_forward, counts.ref_reverse, counts.alt_forward, counts.alt_reverse) == (2, 2, 1, 1)


def test_count_strands_for_indel_uses_mpileup_indel_events():
    counts = count_strands_for_variant(ref="A", alt="AT", read_bases=".+1T,+1t..,")

    assert counts is not None
    assert (counts.ref_forward, counts.ref_reverse, counts.alt_forward, counts.alt_reverse) == (2, 1, 1, 1)


def test_calculate_strand_bias_returns_none_for_symbolic_alleles():
    assert calculate_strand_bias(ref="N", alt="<DEL>", read_bases="") is None


def test_format_pvalue_and_phred_sets_phred_zero_when_pvalue_rounds_to_one():
    formatted_pvalue, formatted_phred = _format_pvalue_and_phred(0.9999999, 7.71462e-15)

    assert formatted_pvalue == "1"
    assert formatted_phred == "0"


def test_calculate_strand_bias_exposes_alt_forward_fraction():
    result = calculate_strand_bias(ref="A", alt="G", read_bases="..,,GGg")

    assert result is not None
    assert result.alt_forward_fraction == 2 / 3


def test_annotate_vcf_strand_bias_writes_info_and_filter(monkeypatch, tmp_path):
    vcf = tmp_path / "input.vcf"
    bam = tmp_path / "sample.bam"
    reference = tmp_path / "reference.fa"
    output = tmp_path / "annotated.vcf"

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tC\tT\t50\tPASS\t.\n",
        encoding="utf-8",
    )
    bam.write_text("bam", encoding="utf-8")
    reference.write_text(
        ">chr1\nC\n",
        encoding="utf-8",
    )

    def fake_mpileup(alignment: Path, reference_path: Path, positions_bed: Path) -> str:
        assert alignment == bam.absolute()
        assert reference_path == reference.absolute()
        bed_lines = positions_bed.read_text(encoding="utf-8").splitlines()
        assert bed_lines == ["chr1\t9\t10"]
        return "chr1\t10\tC\t14\t......,,,,TTTT\t~~~~~~~~~~~~~~\n"

    monkeypatch.setattr("snippy_ng.utils.strand_bias.samtools_mpileup", fake_mpileup)

    annotate_vcf_strand_bias(
        vcf=vcf.absolute(),
        bam=bam.absolute(),
        reference=reference.absolute(),
        output=output,
        filter_pvalue=1.0,
    )

    text = output.read_text(encoding="utf-8")
    assert '##INFO=<ID=SB_REF_FWD' in text
    assert '##INFO=<ID=SB_ALT_FWD_FRAC' in text
    assert '##FORMAT=<ID=SB,' not in text
    assert '##FILTER=<ID=StrandBias' in text
    assert '\tStrandBias\tSB_REF_FWD=6;SB_REF_REV=4;SB_ALT_FWD=4;SB_ALT_REV=0;' in text
    assert 'SB_ALT_FWD_FRAC=1' in text
    assert 'SB_METHOD=FisherExact' in text
