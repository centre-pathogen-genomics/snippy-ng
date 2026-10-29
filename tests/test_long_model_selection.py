import pytest

from snippy_ng.exceptions import Clair3ModelSelectionError
from snippy_ng.pipelines.long import LongPipelineBuilder
from snippy_ng.stages.calling import (
    Clair3Caller,
    LongbowClair3ModelSelector,
    PostDownloadClair3ModelSelector,
)
from snippy_ng.stages.clean_reads import SeqkitCleanLongReads
from snippy_ng.stages.download import DownloadUrlReads
from tests.cli.helpers import make_prepared_reference, stub_load_or_prepare_reference


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


def _build_pipeline(monkeypatch, tmp_path, reads_text):
    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.long.load_or_prepare_reference",
    )
    reads = tmp_path / "reads.fastq"
    reads.write_text(reads_text)
    return LongPipelineBuilder(
        reference=ref_file,
        reads=reads,
        caller="clair3",
        report=False,
        prefix="sample",
    ).build()


def test_build_uses_basecall_header_model_without_longbow(monkeypatch, tmp_path):
    models_root = tmp_path / "models"
    model_dir = models_root / "r1041_e82_400bps_hac_v520"
    model_dir.mkdir(parents=True)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    pipeline = _build_pipeline(
        monkeypatch,
        tmp_path,
        "@read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0\nACGT\n+\nIIII\n",
    )

    assert not any(isinstance(stage, LongbowClair3ModelSelector) for stage in pipeline.stages)
    caller = next(stage for stage in pipeline.stages if isinstance(stage, Clair3Caller))
    assert caller.clair3_model == model_dir.resolve()


def test_build_falls_back_to_longbow_without_basecall_header(monkeypatch, tmp_path):
    models_root = tmp_path / "models"
    (models_root / "r1041_e82_400bps_sup_v520").mkdir(parents=True)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    pipeline = _build_pipeline(monkeypatch, tmp_path, ">read1\nACGT\n")

    assert any(isinstance(stage, LongbowClair3ModelSelector) for stage in pipeline.stages)


def test_build_fails_early_when_header_model_is_missing(monkeypatch, tmp_path):
    models_root = tmp_path / "models"
    (models_root / "r1041_e82_400bps_sup_v520").mkdir(parents=True)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    with pytest.raises(Clair3ModelSelectionError, match="No Clair3 model found"):
        _build_pipeline(
            monkeypatch,
            tmp_path,
            "@read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0\nACGT\n+\nIIII\n",
        )


def test_build_cleans_all_downloaded_read_files(monkeypatch, tmp_path):
    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.long.load_or_prepare_reference",
    )
    r1 = tmp_path / "SRR1_1.fastq.gz"
    r2 = tmp_path / "SRR1_2.fastq.gz"
    r1.touch()
    r2.touch()

    monkeypatch.setattr(
        "snippy_ng.pipelines.long.download_reads",
        lambda accession, stages, **kwargs: [r1, r2],
    )

    pipeline = LongPipelineBuilder(
        reference=ref_file,
        reads_accession="SRR0000001",
        caller="freebayes",
        clean_reads=True,
        report=False,
        prefix="sample",
    ).build()

    cleaner = next(
        stage for stage in pipeline.stages if isinstance(stage, SeqkitCleanLongReads)
    )
    assert cleaner.reads == [r1, r2]


def test_build_feeds_combined_sra_output_to_cleaner(monkeypatch, tmp_path):
    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.long.load_or_prepare_reference",
    )

    pipeline = LongPipelineBuilder(
        reference=ref_file,
        reads_accession="SRR0000001",
        caller="freebayes",
        clean_reads=True,
        report=False,
        prefix="sample",
    ).build()

    cleaner = next(
        stage for stage in pipeline.stages if isinstance(stage, SeqkitCleanLongReads)
    )
    assert len(cleaner.reads) == 1
    assert cleaner.reads[0].name == "SRR0000001.fastq.gz"


def test_build_downloads_url_reads_and_selects_model_after_download(
    monkeypatch, tmp_path
):
    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.long.load_or_prepare_reference",
    )
    models_root = tmp_path / "models"
    (models_root / "r1041_e82_400bps_sup_v520").mkdir(parents=True)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    pipeline = LongPipelineBuilder(
        reference=ref_file,
        reads_url="https://example.com/reads.fastq.gz",
        caller="clair3",
        clean_reads=True,
        report=False,
        prefix="sample",
    ).build()

    download = next(
        stage for stage in pipeline.stages if isinstance(stage, DownloadUrlReads)
    )
    cleaner = next(
        stage for stage in pipeline.stages if isinstance(stage, SeqkitCleanLongReads)
    )
    selectors = [
        stage
        for stage in pipeline.stages
        if isinstance(stage, LongbowClair3ModelSelector)
    ]

    assert cleaner.reads == [download.output.r1]
    assert len(selectors) == 1
    assert type(selectors[0]) is PostDownloadClair3ModelSelector


def test_build_defers_model_selection_until_after_sra_download(monkeypatch, tmp_path):
    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(
        monkeypatch,
        ref_file,
        target="snippy_ng.pipelines.long.load_or_prepare_reference",
    )
    reads = tmp_path / "SRR1_1.fastq.gz"
    reads.touch()
    download_kwargs = {}

    def fake_download_reads(accession, stages, **kwargs):
        download_kwargs.update(kwargs)
        return [reads]

    monkeypatch.setattr(
        "snippy_ng.pipelines.long.download_reads",
        fake_download_reads,
    )
    models_root = tmp_path / "models"
    (models_root / "r1041_e82_400bps_sup_v520").mkdir(parents=True)
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    pipeline = LongPipelineBuilder(
        reference=ref_file,
        reads_accession="SRR0000001",
        caller="clair3",
        report=False,
        prefix="sample",
    ).build()

    selectors = [
        stage
        for stage in pipeline.stages
        if isinstance(stage, LongbowClair3ModelSelector)
    ]
    assert len(selectors) == 1
    assert type(selectors[0]) is PostDownloadClair3ModelSelector
    assert download_kwargs["combine_pairs"] is True


def test_build_fails_early_when_no_models_are_available(monkeypatch, tmp_path):
    models_root = tmp_path / "models"
    models_root.mkdir()
    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    with pytest.raises(Clair3ModelSelectionError, match="No Clair3 models found"):
        _build_pipeline(monkeypatch, tmp_path, ">read1\nACGT\n")
