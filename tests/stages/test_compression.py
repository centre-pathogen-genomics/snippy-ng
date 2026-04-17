from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.compression import CramCompressor


def test_cram_compressor_does_not_redirect_stdout_to_output():
    stage = CramCompressor(
        input=Path("sample.filtered.bam"),
        reference=Path("reference.fa"),
    )

    commands = stage.create_commands(Context(cpus=4))

    assert len(commands) == 1
    command = commands[0]
    assert command.output_file is None
    assert command.command == [
        "samtools",
        "view",
        "--threads", "4",
        "-C",
        "-T",
        "reference.fa",
        "--output-fmt-option",
        "embed_ref=1",
        "-o",
        "sample.filtered.cram",
        "sample.filtered.bam",
    ]
