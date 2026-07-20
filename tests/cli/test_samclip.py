from click.testing import CliRunner

from snippy_ng.cli.samclip_cli import samclip


def test_samclip_fraction_option_replaces_default_absolute_limit(tmp_path):
    index = tmp_path / "ref.fa.fai"
    sam = tmp_path / "input.sam"
    index.write_text("ref\t1000\t0\t80\t81\n")
    sam.write_text(
        "@HD\tVN:1.6\tSO:queryname\n"
        "@SQ\tSN:ref\tLN:1000\n"
        f"r1\t0\tref\t500\t60\t50S50M\t*\t0\t0\t{'A' * 100}\t{'I' * 100}\n"
    )

    result = CliRunner().invoke(
        samclip,
        ["--index", str(index), "--max-clip-fraction", "0.5", str(sam)],
    )

    assert result.exit_code == 0, result.output
    assert "r1\t0\tref" in result.output
