from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.alignment import AssemblyAligner, AssemblyNucmerAligner, Minimap2ShortReadAligner


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
    assert commands[2][1:5] == ["-m", "snippy_ng", "utils", "aln"]
    assert commands[2][5] == "samclip"
    assert "--fix-mate" in commands[2]
    assert commands[3][:2] == ["samtools", "fixmate"]
    assert commands[4][:2] == ["samtools", "sort"]
    assert commands[5][:2] == ["samtools", "markdup"]


def test_nucmer_assembly_aligner_honours_configured_tunables(tmp_path):
    stage = AssemblyNucmerAligner(
        reference=Path("reference.fa"),
        assembly=Path("assembly.fa"),
        breaklen=250,
        mincluster=80,
        maxgap=120,
        minmatch=24,
        minalign=500,
        prefix="sample",
    )

    command = stage.create_commands(Context(cpus=4, ram=8, tmpdir=tmp_path))[0].command

    assert command == [
        "nucmer",
        "--prefix",
        "sample",
        "--threads",
        "4",
        "--breaklen",
        "250",
        "--mincluster",
        "80",
        "--maxgap",
        "120",
        "--minmatch",
        "24",
        "--minalign",
        "500",
        "reference.fa",
        "assembly.fa",
    ]


def test_minimap2_assembly_aligner_uses_configured_preset(tmp_path):
    stage = AssemblyAligner(
        reference=Path("reference.fa"),
        assembly=Path("assembly.fa"),
        minimap_preset="asm5",
        prefix="sample",
    )

    pipeline = stage.create_commands(Context(cpus=4, ram=8, tmpdir=tmp_path))[0]
    commands = [command.command for command in pipeline.processes]

    assert commands[0] == [
        "minimap2",
        "-x",
        "asm5",
        "-t",
        "4",
        "-c",
        "--cs",
        "reference.fa",
        "assembly.fa",
    ]
    assert commands[1] == ["sort", "-k6,6", "-k8,8n"]
