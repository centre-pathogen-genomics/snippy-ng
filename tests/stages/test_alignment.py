from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.alignment import Minimap2ShortReadAligner


def test_minimap2_short_read_pipeline_name_sorts_before_filtering(tmp_path):
    stage = Minimap2ShortReadAligner(
        reference=Path("reference.fa"),
        reads=[Path("reads_1.fq.gz"), Path("reads_2.fq.gz")],
        prefix="sample",
    )

    pipeline = stage.create_commands(Context(cpus=4, ram=8, tmpdir=tmp_path))[0]
    commands = [command.command for command in pipeline.processes]

    assert commands[0][:4] == ["minimap2", "-a", "-x", "sr"]
    assert commands[1][:5] == ["samtools", "sort", "-n", "-O", "sam"]
    assert commands[2][1:4] == ["-m", "snippy_ng", "utils"]
    assert commands[2][4] == "samclip"
    assert "--fix-mate" in commands[2]
    assert commands[3][:2] == ["samtools", "fixmate"]
    assert commands[4][:2] == ["samtools", "sort"]
    assert commands[5][:2] == ["samtools", "markdup"]
