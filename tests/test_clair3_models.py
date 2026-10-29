import gzip

import pytest

from snippy_ng.exceptions import Clair3ModelSelectionError
from snippy_ng.utils.clair3_models import choose_model, get_basecall_model_id


@pytest.fixture(autouse=True)
def isolated_model_roots(monkeypatch):
    for env_var in (
        "CLAIR3_MODELS",
        "CLAIR3_MODEL_ROOT",
        "CLAIR3_MODEL_CONTAINER_ROOT",
        "CLAIR3_MODELS_CONTAINER_ROOT",
        "CONDA_PREFIX",
    ):
        monkeypatch.delenv(env_var, raising=False)


def _write_reads(path, header="read1", read_id="dna_r10.4.1_e8.2_400bps_hac@v5.2.0"):
    path.write_text(f"@{header} basecall_model_version_id={read_id}\nACGT\n+\nIIII\n")


def _make_model(root, name):
    model_dir = root / name
    model_dir.mkdir(parents=True)
    return model_dir


def test_choose_model_from_key_value_header(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected = _make_model(models_root, "r1041_e82_400bps_hac_v520")
    _make_model(models_root, "r1041_e82_400bps_sup_v520")
    reads = tmp_path / "reads.fastq"
    _write_reads(reads)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    assert choose_model(reads) == expected.resolve()


def test_choose_model_from_gzipped_sam_tag_header(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected = _make_model(models_root, "r1041_e82_400bps_sup_v500")
    reads = tmp_path / "reads.fastq.gz"
    with gzip.open(reads, "wt") as handle:
        handle.write(
            "@read1 RG:Z:1234_dna_r10.4.1_e8.2_400bps_sup@v5.0.0_barcode01\nACGT\n+\nIIII\n"
        )
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    assert choose_model(reads) == expected.resolve()


def test_choose_model_returns_none_without_basecall_id(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    _make_model(models_root, "r1041_e82_400bps_sup_v520")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@read1\nACGT\n+\nIIII\n")
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    assert choose_model(reads) is None


def test_choose_model_errors_when_models_dir_has_no_models(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    models_root.mkdir()
    reads = tmp_path / "reads.fastq"
    _write_reads(reads)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    with pytest.raises(Clair3ModelSelectionError, match="No Clair3 models found"):
        choose_model(reads)


def test_choose_model_errors_when_model_is_missing(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    _make_model(models_root, "r1041_e82_400bps_sup_v520")
    reads = tmp_path / "reads.fastq"
    _write_reads(reads)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    with pytest.raises(Clair3ModelSelectionError, match="No Clair3 model found") as exc_info:
        choose_model(reads)

    assert str(models_root.resolve()) in str(exc_info.value)
    assert "No model matches preset 'hac'" in str(exc_info.value)
    assert (
        "snippy-ng utils model download r1041_e82_400bps_hac_v520"
        in str(exc_info.value)
    )


def test_get_basecall_model_id_supports_both_header_formats():
    assert (
        get_basecall_model_id(
            "read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
        )
        == "dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
    )
    assert (
        get_basecall_model_id(
            "read1 RG:Z:uuid_dna_r10.4.1_e8.2_400bps_hac@v5.2.0_barcode01"
        )
        == "dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
    )
    assert get_basecall_model_id("read1") is None


def test_basecall_model_id_to_model_name():
    from snippy_ng.utils.clair3_models import basecall_model_id_to_model_name

    assert (
        basecall_model_id_to_model_name("dna_r10.4.1_e8.2_400bps_hac@v4.2.0")
        == "r1041_e82_400bps_hac_v420"
    )
    assert (
        basecall_model_id_to_model_name("r1041_e82_400bps_sup_v520")
        == "r1041_e82_400bps_sup_v520"
    )


def test_download_clair3_model_writes_files(tmp_path, monkeypatch):
    import io
    import urllib.request

    from snippy_ng.utils.clair3_models import download_clair3_model

    def fake_urlopen(request: urllib.request.Request):
        return io.BytesIO(f"model:{request.full_url}".encode())

    monkeypatch.setattr("snippy_ng.utils.clair3_models.urllib.request.urlopen", fake_urlopen)

    model_dir = download_clair3_model(
        "dna_r10.4.1_e8.2_400bps_hac@v4.2.0", tmp_path
    )

    assert model_dir == tmp_path / "r1041_e82_400bps_hac_v420"
    assert (model_dir / "pileup.pt").read_bytes() == (
        b"model:https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/"
        b"r1041_e82_400bps_hac_v420/pileup.pt"
    )
    assert (model_dir / "full_alignment.pt").is_file()


def test_download_clair3_model_skips_existing_files(tmp_path, monkeypatch):
    from snippy_ng.utils.clair3_models import download_clair3_model

    model_dir = tmp_path / "r1041_e82_400bps_sup_v520"
    model_dir.mkdir(parents=True)
    (model_dir / "pileup.pt").write_bytes(b"pileup")
    (model_dir / "full_alignment.pt").write_bytes(b"full")

    def fail_urlopen(request):
        raise AssertionError("download should be skipped")

    monkeypatch.setattr("snippy_ng.utils.clair3_models.urllib.request.urlopen", fail_urlopen)

    assert download_clair3_model("r1041_e82_400bps_sup_v520", tmp_path) == model_dir


def test_download_clair3_model_reports_missing_model(tmp_path, monkeypatch):
    import urllib.error

    from snippy_ng.exceptions import Clair3ModelDownloadError
    from snippy_ng.utils.clair3_models import download_clair3_model

    def missing_urlopen(request):
        raise urllib.error.HTTPError(request.full_url, 404, "Not Found", None, None)

    monkeypatch.setattr("snippy_ng.utils.clair3_models.urllib.request.urlopen", missing_urlopen)

    with pytest.raises(Clair3ModelDownloadError, match="is not available"):
        download_clair3_model("dna_r10.4.1_e8.2_400bps_hac@v4.2.0", tmp_path)


def test_basecall_model_id_to_model_name_rejects_path_traversal():
    from snippy_ng.exceptions import Clair3ModelDownloadError
    from snippy_ng.utils.clair3_models import basecall_model_id_to_model_name

    with pytest.raises(Clair3ModelDownloadError, match="not a valid"):
        basecall_model_id_to_model_name("..")
