from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.calling import (
    FreebayesCaller,
    MIN_FREEBAYES_CHUNK_SIZE,
    get_freebayes_chunk_size,
)


def test_get_freebayes_chunk_size_scales_with_cpus_and_reference_index(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000000\t0\t0\t0\n")

    num_chunks, chunk_size = get_freebayes_chunk_size(reference, reference_index, cpus=4)

    assert num_chunks == 7
    assert chunk_size == 142857


def test_get_freebayes_chunk_size_honours_minimum(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000\t0\t0\t0\n")

    _, chunk_size = get_freebayes_chunk_size(reference, reference_index, cpus=8)

    assert chunk_size == MIN_FREEBAYES_CHUNK_SIZE


def test_freebayes_caller_uses_adaptive_chunk_size_for_region_generation(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000000\t0\t0\t0\n")
    bam = tmp_path / "reads.bam"
    bam.write_text("")
    bam_index = tmp_path / "reads.bam.bai"
    bam_index.write_text("")

    stage = FreebayesCaller(
        reference=reference,
        reference_index=reference_index,
        bam=bam,
        bam_index=bam_index,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=4))

    generate_regions = commands[0].processes[0]
    assert generate_regions.command == [
        "fasta_generate_regions.py",
        str(reference_index),
        "142857",
    ]
