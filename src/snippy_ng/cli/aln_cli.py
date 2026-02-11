import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, GlobalOption, add_snippy_global_options, create_outdir_callback
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=Path("core"), required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=Path), help="Output directory for the prepared reference", callback=create_outdir_callback, cls=GlobalOption)
@add_snippy_global_options(exclude=['prefix', 'outdir'])
@click.argument("snippy_dirs", required=True, nargs=-1, type=click.Path(exists=True, resolve_path=True, readable=True))
@click.option("--ref", "reference", type=click.Path(exists=True, resolve_path=True, readable=True), required=True, help="Reference FASTA used to define contig order")
@click.option("--core", type=click.FLOAT, default=0.95, help="Proportion of samples a site must be present in to be included in the core alignment (0.0-1.0)")
def aln(**config):
    """
    Create alignment from multiple snippy runs
    """
    from snippy_ng.pipelines.aln import create_aln_pipeline_stages
    from snippy_ng.pipelines.pipeline_runner import run_snippy_pipeline

    if not config.get("snippy_dirs"):
        raise click.UsageError("Please provide at least one snippy directory!")

    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    stages = create_aln_pipeline_stages(
        snippy_dirs=config["snippy_dirs"],
        reference=Path(config["reference"]),
        core=config["core"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    )

    # Run the pipeline
    return run_snippy_pipeline(
        stages,
        skip_check=config['skip_check'],
        check=config['check'],
        outdir=config['outdir'],
        quiet=config['quiet'],
        create_missing=config['create_missing'],
        keep_incomplete=config['keep_incomplete'],
    )