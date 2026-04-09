from pathlib import Path
import os

import pytest

from snippy_ng.context import Context
from snippy_ng.stages.calling import (
    Clair3Caller,
    Clair3ModelSelectorError,
    FreebayesCaller,
    LongbowClair3ModelSelector,
    MIN_FREEBAYES_CHUNK_SIZE,
    get_short_chunk_size,
)


def test_get_short_chunk_size_scales_with_cpus_and_reference_index(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000000\t0\t0\t0\n")

    num_chunks, chunk_size = get_short_chunk_size(reference, reference_index, cpus=4)

    assert num_chunks == 7
    assert chunk_size == 142857


def test_get_short_chunk_size_honours_minimum(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000\t0\t0\t0\n")

    _, chunk_size = get_short_chunk_size(reference, reference_index, cpus=8)

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


def test_longbow_resolver_picks_r10_sup_model_from_env_root(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected_model = models_root / "r1041_e82_400bps_sup_v520"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell_version":"R10","basecaller":"Dorado","major_version":"Dorado0","basecalling_mode":"SUP","dorado_model_version":"V5.0.0"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_uses_launch_dir_for_relative_env_root(tmp_path, monkeypatch):
    launch_dir = tmp_path / "launch"
    launch_dir.mkdir()
    models_root = launch_dir / "clair3_models"
    expected_model = models_root / "r1041_e82_400bps_sup_v520"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"Sample":"reads.fastq.gz","Flowcell":"R10","Software":"guppy","Version":"5or6","Mode":"SUP"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("PWD", str(launch_dir))
    monkeypatch.setenv("CLAIR3_MODELS", "./clair3_models")

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_picks_r9_guppy5_legacy_model(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected_model = models_root / "ont_guppy5"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell":"R9","basecaller":"Guppy","major_version":"Guppy5/6","mode":"HAC"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_parses_title_case_longbow_json(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected_model = models_root / "r1041_e82_400bps_sup_v520"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"Sample":"reads.fastq.gz","Flowcell":"R10","Software":"guppy","Version":"5or6","Mode":"SUP","Confidence level":"very high"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_errors_when_no_model_root_exists(tmp_path, monkeypatch):
    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell":"R10","basecaller":"Dorado","major_version":"Dorado0","mode":"SUP"}'
    )

    monkeypatch.delenv("CLAIR3_MODELS", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_ROOT", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_PATH", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_CONTAINER_ROOT", raising=False)
    monkeypatch.delenv("CLAIR3_MODELS_CONTAINER_ROOT", raising=False)
    monkeypatch.delenv("CONDA_PREFIX", raising=False)

    with pytest.raises(Clair3ModelSelectorError, match="Could not find any Clair3 model roots"):
        LongbowClair3ModelSelector.resolve_clair3_model(
            prediction_json,
            tmp_path / "resolved_model",
        )


def test_longbow_resolver_allows_container_only_model_root(tmp_path, monkeypatch):
    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell_version":"R10","basecaller":"Dorado","major_version":"Dorado0","basecalling_mode":"SUP","dorado_model_version":"V5.0.0"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.delenv("CLAIR3_MODELS", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_ROOT", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_PATH", raising=False)
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    monkeypatch.setenv("CLAIR3_MODEL_CONTAINER_ROOT", "/opt/models")

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == "/opt/models/r1041_e82_400bps_sup_v500"


def test_clair3_caller_reads_model_path_from_manifest_file(tmp_path):
    manifest = tmp_path / "snippy.clair3_model.txt"
    manifest.write_text("/opt/models/r1041_e82_400bps_sup_v520\n", encoding="utf-8")
    bam = tmp_path / "reads.bam"
    bam.write_text("")
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1\t0\t1\t2\n")

    stage = Clair3Caller(
        bam=bam,
        reference=reference,
        reference_index=reference_index,
        clair3_model=manifest,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=1))

    assert commands[0].command[1] == "--model_path=/opt/models/r1041_e82_400bps_sup_v520"
    assert "--chunk_size=10000" in commands[0].command
