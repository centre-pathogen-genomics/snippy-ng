import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import click

from snippy_ng.cli.utils import (
    is_sra_accession,
    is_url,
    reads_or_url_or_accession_callback,
)
from snippy_ng.exceptions import Clair3ModelDownloadError, SnippyError


def default_model_directory() -> Path:
    """Return the configured Clair3 model root, or ./clair3_models."""
    for env_var in ("CLAIR3_MODELS", "CLAIR3_MODEL_ROOT"):
        value = os.environ.get(env_var)
        if value:
            return Path(value).expanduser()
    return Path("clair3_models")


@click.group()
def model():
    """Clair3 model utilities."""


@model.command("download")
@click.argument("model", metavar="MODEL")
@click.option(
    "--outdir",
    "-o",
    default=None,
    type=click.Path(file_okay=False, path_type=Path),
    help="Model directory root (default: $CLAIR3_MODELS or ./clair3_models)",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-download the model even if it already exists",
)
def download(model: str, outdir: Optional[Path], force: bool):
    """
    Download a Clair3 model.

    MODEL is a basecall model id or a Clair3 model directory name, for example
    dna_r10.4.1_e8.2_400bps_sup@v5.2.0 or r1041_e82_400bps_sup_v520.
    """
    from snippy_ng.utils.clair3_models import download_clair3_model

    destination = Path(outdir).absolute() if outdir is not None else default_model_directory()
    try:
        model_dir = download_clair3_model(model, destination, force=force)
    except Clair3ModelDownloadError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Clair3 model available at: {model_dir}")


@model.command("check")
@click.argument(
    "reads",
    metavar="READS",
    type=click.STRING,
    callback=reads_or_url_or_accession_callback,
)
@click.option(
    "--outdir",
    "-o",
    default=Path("."),
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory for Longbow prediction and model selection outputs",
)
@click.option(
    "--prefix",
    "-p",
    default="snippy",
    help="Prefix for Longbow prediction and model selection outputs",
)
@click.option(
    "--cpus",
    default=1,
    type=click.IntRange(min=1),
    help="Threads passed to Longbow when it is used",
)
def check(reads: str, outdir: Path, prefix: str, cpus: int):
    """
    Select the Clair3 model for READS.

    The read header is checked first. If it has no basecall model id, Longbow is
    run and the predicted model is resolved from the configured model roots.
    """
    from snippy_ng.stages.calling import PostDownloadClair3ModelSelector
    from snippy_ng.stages.download import DownloadUrlReads

    if is_sra_accession(reads):
        raise click.ClickException(
            "Model check requires a local reads file or an HTTP(S) URL, not an SRA accession."
        )

    outdir = Path(outdir).absolute()
    outdir.mkdir(parents=True, exist_ok=True)
    prediction_json = outdir / f"{prefix}.longbow.json"
    model_file = outdir / f"{prefix}.clair3_model.txt"

    try:
        if is_url(reads):
            with tempfile.TemporaryDirectory(prefix="snippy-ng-model-check-") as tmpdir:
                download_stage = DownloadUrlReads(
                    url=reads,
                    output_directory=Path(tmpdir),
                )
                download_stage.download_file(reads, download_stage.output.r1)
                PostDownloadClair3ModelSelector.select_from_header_or_predict(
                    download_stage.output.r1,
                    prediction_json,
                    model_file,
                    cpus,
                )
        else:
            PostDownloadClair3ModelSelector.select_from_header_or_predict(
                Path(reads),
                prediction_json,
                model_file,
                cpus,
            )
        model_dir = model_file.read_text(encoding="utf-8").strip()
    except SnippyError as exc:
        raise click.ClickException(str(exc)) from exc
    except subprocess.CalledProcessError as exc:
        raise click.ClickException(
            f"Longbow failed with exit code {exc.returncode}."
        ) from exc

    click.echo(f"Clair3 model: {model_dir}")
