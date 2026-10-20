import click
from typing import Any, Optional
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import (
    CommandWithGlobals,
    GlobalOption,
    add_snippy_global_options,
    check_outdir_callback,
)
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={"show_default": True})
@click.option(
    "--outdir",
    "-o",
    default=Path("core"),
    required=False,
    type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True),
    help="Output directory for the prepared reference",
    callback=check_outdir_callback,
    cls=GlobalOption,
)
@click.option(
    "--prefix", "-p", default="core", help="Prefix for output files", cls=GlobalOption
)
@add_snippy_global_options(exclude=["prefix", "outdir", "sample_name"])
@click.argument(
    "snippy_dirs",
    required=True,
    nargs=-1,
    type=AbsolutePath(exists=True, readable=True),
)
@click.option(
    "--reference",
    "--ref",
    type=AbsolutePath(exists=True, readable=True),
    required=True,
    help="Reference FASTA used to define contig order",
)
@click.option(
    "--core",
    type=click.FloatRange(min=0, max=1.0),
    default=0.95,
    help="Proportion of samples a site must be present in to be included in the core alignment",
)
@click.option(
    "--remove-alignment-outliers",
    "remove_outliers",
    is_flag=True,
    default=False,
    help="Remove detected low-alignment outlier samples from the alignment before core filtering",
)
@click.option(
    "--alignment-outlier-threshold",
    "outlier_threshold",
    type=click.FloatRange(min=0.0, max=1.0, min_open=True),
    default=0.20,
    help="Posterior probability threshold for retaining membership in the retained alignment percentage clusters. Only used with --remove-alignment-outliers and cohorts large enough for mixture modelling",
)
@click.option(
    "--snp-distance-format",
    type=click.Choice(["tabular", "phylip"]),
    default="tabular",
    help="Output format for pairwise snip distance matrix",
)
@click.option(
    "--mask",
    default=None,
    type=AbsolutePath(exists=True, readable=True),
    help="Mask file (BED format, reference coordinates) to mask regions with Ns in every sample of the alignment",
)
def core(
    snippy_dirs: tuple[Path, ...],
    reference: Path,
    core: float,
    remove_outliers: bool,
    outlier_threshold: float,
    snp_distance_format: str,
    mask: Optional[Path],
    outdir: Path,
    prefix: str,
    **context: Any,
):
    """
    Create core alignment from multiple Snippy-NG runs
    """
    from snippy_ng.context import Context
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
        remove_outliers=remove_outliers,
        outlier_threshold=outlier_threshold,
        snp_distance_format=snp_distance_format,
        mask=mask,
        prefix=prefix,
    ).build()

    # Run the pipeline
    context["outdir"] = outdir
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
