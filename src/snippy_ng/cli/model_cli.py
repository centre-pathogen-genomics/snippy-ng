import os
from pathlib import Path
from typing import Optional

import click

from snippy_ng.exceptions import Clair3ModelDownloadError


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
