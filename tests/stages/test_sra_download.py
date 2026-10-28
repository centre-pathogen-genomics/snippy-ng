import io
from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages import PythonCommand
from snippy_ng.stages.download import DownloadSraReads, DownloadUrlReads


def test_download_url_reads_writes_file(tmp_path: Path, monkeypatch):
    stage = DownloadUrlReads(
        url="https://example.com/data/sample.fastq.gz?token=abc",
        output_directory=tmp_path,
    )
    monkeypatch.setattr(
        "snippy_ng.stages.download.urllib.request.urlopen",
        lambda request: io.BytesIO(b"@read1\nACGT\n+\nIIII\n"),
    )

    stage.download_file(stage.url, stage.output.r1)

    assert stage.output.r1 == tmp_path / "sample.fastq.gz"
    assert stage.output.r1.read_bytes() == b"@read1\nACGT\n+\nIIII\n"


def test_download_sra_reads_combines_downloaded_files(tmp_path: Path):
    stage = DownloadSraReads(
        accession="SRR123456",
        output_directory=tmp_path,
        combine_pairs=True,
    )
    (tmp_path / "SRR123456_1.fastq.gz").write_bytes(b"r1\n")
    (tmp_path / "SRR123456_2.fastq.gz").write_bytes(b"r2\n")

    stage._combine_sracha_output()

    assert stage.output.r1 == tmp_path / "SRR123456.fastq.gz"
    assert stage.output.r2 is None
    assert stage.output.r1.read_bytes() == b"r1\nr2\n"


def test_sracha_normalize_python_command_accepts_no_context(tmp_path: Path):
    stage = DownloadSraReads(
        accession="SRR123456",
        output_directory=tmp_path,
    )
    single_end = tmp_path / "SRR123456.fastq.gz"
    r1 = tmp_path / "SRR123456_1.fastq.gz"
    single_end.write_bytes(b"reads")

    command = stage.create_commands(Context())[1]

    assert isinstance(command, PythonCommand)
    command.func(*command.args)

    assert not single_end.exists()
    assert r1.read_bytes() == b"reads"

