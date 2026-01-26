import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options, create_outdir_callback
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@snippy_global_options
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--outdir", "-o", default=Path("reference"), required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=Path), help="Output directory for the prepared reference", callback=create_outdir_callback)
def ref(**config):
    """
    Reference preparation pipeline

    Examples:

        $ snippy-ng ref --reference ref.fa --outdir output
    """
    from snippy_ng.pipelines.common import prepare_reference
    from snippy_ng.cli.utils.pipeline_runner import run_snippy_pipeline

    run_snippy_pipeline(
        stages=[
            prepare_reference(
                reference_path=config["reference"],
                output_directory=config["outdir"]
            )
        ],
        skip_check=config["skip_check"],
        check=config["check"],
        outdir=Path("."),
        quiet=config["quiet"],
        continue_last_run=config["continue_last_run"],
        keep_incomplete=config["keep_incomplete"],
    )