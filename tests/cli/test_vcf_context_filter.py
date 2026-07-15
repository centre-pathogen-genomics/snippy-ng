from click.testing import CliRunner

from snippy_ng.cli.vcf.context_filter_cli import context_filter


def test_vcf_context_filter_cli_marks_nearby_snp(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t10\t.\tA\tT\t60\t.\t.",
                "chr1\t15\t.\tAC\tA\t60\t.\t.",
            ]
        ) + "\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        context_filter,
        [
            str(input_vcf),
            "--output",
            str(output_vcf),
            "--min-snp-distance-to-indel",
            "10",
        ],
    )

    assert result.exit_code == 0, result.output
    records = [line for line in output_vcf.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    assert records[0].split("\t")[6] == "LowQual"
    assert "CONTEXT_LOWQUAL_REASON=NEAR_INDEL" in records[0]
