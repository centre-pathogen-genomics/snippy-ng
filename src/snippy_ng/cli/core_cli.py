import click
from typing import Any
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, GlobalOption, add_snippy_global_options, create_outdir_callback
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=Path("core"), required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True), help="Output directory for the prepared reference", callback=create_outdir_callback, cls=GlobalOption)
@click.option("--prefix", "-p", default="core", help="Prefix for output files", cls=GlobalOption)
@add_snippy_global_options(exclude=['prefix', 'outdir'])
@click.argument("snippy_dirs", required=True, nargs=-1, type=AbsolutePath(exists=True, readable=True))
@click.option("--ref", "reference", type=AbsolutePath(exists=True, readable=True), required=True, help="Reference FASTA used to define contig order")
@click.option("--core", type=click.FLOAT, default=0.95, help="Proportion of samples a site must be present in to be included in the core alignment (0.0-1.0)")
def core(snippy_dirs: tuple[Path, ...], reference: Path, core: float, outdir: Path, prefix: str, **context: Any):
    """
    Create core alignment from multiple Snippy-NG runs
    """
    from snippy_ng import Context
    from snippy_ng.pipelines.core import CorePipelineBuilder

    if not snippy_dirs:
        raise click.UsageError("Please provide at least one snippy directory!")

    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = CorePipelineBuilder(
        snippy_dirs=snippy_dirs,
        reference=Path(reference),
        core=core,
        prefix=prefix,
    ).build()

    # Run the pipeline
    context["outdir"] = outdir
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)