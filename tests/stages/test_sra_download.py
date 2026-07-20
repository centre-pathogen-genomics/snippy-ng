from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages import PythonCommand
from snippy_ng.stages.download import DownloadSraReads


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

