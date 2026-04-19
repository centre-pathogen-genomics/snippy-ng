import subprocess

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng


def test_cnv_cli_outputs_copy_number_table(monkeypatch, tmp_path):
    alignment = tmp_path / "sample.cram"
    alignment.write_text("cram")
    captured = {}

    def fake_run(command, check, stdout, stderr, text):
        captured["command"] = command
        assert check is True
        assert stdout == subprocess.PIPE
        assert stderr == subprocess.PIPE
        assert text is True
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=(
                "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
                "chr1\t1\t1000\t100\t1000\t100\t30\t40\t60\n"
                "plasmid\t1\t100\t100\t100\t100\t90\t40\t60\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.cnv.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "cnv", str(alignment)],
    )

    assert result.exit_code == 0, result.output
    assert captured["command"] == [
        "samtools",
        "coverage",
        str(alignment.resolve()),
    ]
    assert result.output == (
        "contig_id\tread_depth\tcopy_number\n"
        "chr1\t30\t1\n"
        "plasmid\t90\t3\n"
    )


def test_cnv_cli_no_header(monkeypatch, tmp_path):
    alignment = tmp_path / "sample.bam"
    alignment.write_text("bam")

    def fake_run(command, check, stdout, stderr, text):
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=(
                "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
                "chr1\t1\t1000\t100\t1000\t100\t30\t40\t60\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.cnv.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "cnv", "--no-header", str(alignment)],
    )

    assert result.exit_code == 0, result.output
    assert result.output == "chr1\t30\t1\n"


def test_cnv_cli_gff_outputs_feature_copy_number_table(monkeypatch, tmp_path):
    alignment = tmp_path / "sample.cram"
    alignment.write_text("cram")
    gff = tmp_path / "reference.gff"
    gff.write_text(
        "chr1\t.\tCDS\t1\t3\t.\t+\t0\tID=cds1\n"
        "plasmid\t.\tCDS\t1\t3\t.\t+\t0\tID=cds2\n"
    )
    commands = []

    def fake_run(command, check, stdout, stderr, text):
        commands.append(command)
        if command[1] == "coverage":
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=(
                    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
                    "chr1\t1\t1000\t100\t1000\t100\t30\t40\t60\n"
                    "plasmid\t1\t100\t100\t100\t100\t90\t40\t60\n"
                ),
                stderr="",
            )
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=(
                "chr1\t1\t29\n"
                "chr1\t2\t30\n"
                "chr1\t3\t1000\n"
                "plasmid\t1\t89\n"
                "plasmid\t2\t90\n"
                "plasmid\t3\t91\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.cnv.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "cnv",
            str(alignment),
            "--gff",
            str(gff),
            "--feature",
            "CDS",
        ],
    )

    assert result.exit_code == 0, result.output
    assert commands[0] == [
        "samtools",
        "coverage",
        str(alignment.resolve()),
    ]
    assert commands[1][:5] == ["samtools", "depth", "-aa", "-b", commands[1][4]]
    assert commands[1][5:] == [str(alignment.resolve())]
    assert result.output == (
        "feature_id\tcontig_id\tstart\tend\tread_depth\tcopy_number\n"
        "cds1\tchr1\t1\t3\t30\t1\n"
        "cds2\tplasmid\t1\t3\t90\t3\n"
    )


def test_cnv_cli_known_single_copy_overrides_contig_baseline(monkeypatch, tmp_path):
    alignment = tmp_path / "sample.cram"
    alignment.write_text("cram")
    commands = []

    def fake_run(command, check, stdout, stderr, text):
        commands.append(command)
        if command[1] == "coverage":
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=(
                    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
                    "chr1\t1\t1000\t100\t1000\t100\t60\t40\t60\n"
                    "plasmid\t1\t100\t100\t100\t100\t120\t40\t60\n"
                ),
                stderr="",
            )
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=(
                "chr1\t10\t29\n"
                "chr1\t11\t30\n"
                "chr1\t12\t31\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.cnv.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "cnv",
            str(alignment),
            "--known-single-copy",
            "10,12",
        ],
    )

    assert result.exit_code == 0, result.output
    assert commands[1] == [
        "samtools",
        "depth",
        "-aa",
        "-r",
        "chr1:10-12",
        str(alignment.resolve()),
    ]
    assert result.output == (
        "contig_id\tread_depth\tcopy_number\n"
        "chr1\t60\t2\n"
        "plasmid\t120\t4\n"
    )


def test_cnv_cli_known_single_copy_named_region_overrides_feature_baseline(monkeypatch, tmp_path):
    alignment = tmp_path / "sample.cram"
    alignment.write_text("cram")
    gff = tmp_path / "reference.gff"
    gff.write_text("plasmid\t.\tCDS\t1\t3\t.\t+\t0\tID=cds1\n")
    commands = []

    def fake_run(command, check, stdout, stderr, text):
        commands.append(command)
        if command[1] == "coverage":
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=(
                    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
                    "chr1\t1\t1000\t100\t1000\t100\t60\t40\t60\n"
                    "plasmid\t1\t100\t100\t100\t100\t120\t40\t60\n"
                ),
                stderr="",
            )
        if "-r" in command:
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=(
                    "chr1\t10\t29\n"
                    "chr1\t11\t30\n"
                    "chr1\t12\t31\n"
                ),
                stderr="",
            )
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=(
                "plasmid\t1\t119\n"
                "plasmid\t2\t120\n"
                "plasmid\t3\t121\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.cnv.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "cnv",
            str(alignment),
            "--gff",
            str(gff),
            "--feature",
            "CDS",
            "--known-single-copy",
            "chr1:10-12",
        ],
    )

    assert result.exit_code == 0, result.output
    assert commands[1][0:5] == ["samtools", "depth", "-aa", "-r", "chr1:10-12"]
    assert result.output == (
        "feature_id\tcontig_id\tstart\tend\tread_depth\tcopy_number\n"
        "cds1\tplasmid\t1\t3\t120\t4\n"
    )
