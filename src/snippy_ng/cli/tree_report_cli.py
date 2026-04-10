import click
from typing import Any, Optional
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, GlobalOption, add_snippy_global_options, check_outdir_callback
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default="report", required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=Path), help="Output directory for phylogenetic tree results", callback=check_outdir_callback, cls=GlobalOption)
@click.option("--prefix", "-p", default="report", help="Prefix for html report", cls=GlobalOption)
@add_snippy_global_options(exclude=['prefix', 'outdir'])
@click.argument("newick", required=True, type=AbsolutePath(exists=True, readable=True))
@click.option("--metadata", required=False, type=AbsolutePath(exists=True, readable=True), help="Optional metadata file (JSON or CSV) to include in the report")
@click.option("--color-by-column", required=False, type=click.STRING, help="Column name in the metadata to color the tree by")
@click.option("--logs", required=False, type=AbsolutePath(exists=True, readable=True), help="Optional log file to include in the report")
@click.option("--title", required=False, type=click.STRING, default="Snippy-NG Report", help="Title for the HTML report")
@click.option("--mid-point-root", is_flag=True, default=False, help="Mid-point root the tree in the report")
@click.option("--ladderize", is_flag=True, default=False, help="Ladderize the tree in the report")
def tree_report(
    newick: Path,
    metadata: Optional[Path],
    color_by_column: Optional[str],
    logs: Optional[Path],
    title: str,
    outdir: Optional[Path],
    prefix: str,
    mid_point_root: bool,
    ladderize: bool,
    **context: Any,
):
    """
    Create an HTML report with an interactive phylogenetic tree

    Example usage:

        snippy-ng utils report tree.newick --metadata metadata.csv
    """
    from snippy_ng.context import Context
    from snippy_ng.pipelines.report import ReportPipelineBuilder
   
    if metadata and metadata.suffix.lower() not in [".json", ".csv", ".tsv"]:
        raise click.BadParameter("Metadata file must be in JSON, CSV, or TSV format", param_hint="--metadata")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = ReportPipelineBuilder(
        tree=newick,
        mid_point_root=mid_point_root,
        ladderize=ladderize,
        title=title,
        metadata=metadata,
        color_by_column=color_by_column,
        logs=logs,
        prefix=prefix,
    ).build()

    context["outdir"] = outdir or Path(".").absolute()
    # Run the pipeline
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
