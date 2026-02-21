import click
from typing import Any, Optional
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, GlobalOption, add_snippy_global_options, create_outdir_callback
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=None, required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True), help="Output directory for phylogenetic tree results", callback=create_outdir_callback, cls=GlobalOption)
@click.option("--prefix", "-p", default="report", help="Prefix for html report", cls=GlobalOption)
@add_snippy_global_options(exclude=['prefix', 'outdir'])
@click.option("--tree", required=True, type=AbsolutePath(exists=True, readable=True), help="Newick tree file to include in the report")
@click.option("--metadata", required=False, type=AbsolutePath(exists=True, readable=True), help="Optional metadata file (JSON or CSV) to include in the report")
@click.option("--logs", required=False, type=AbsolutePath(exists=True, readable=True), help="Optional log file to include in the report")
@click.option("--title", required=False, type=click.STRING, default="Snippy-NG Report", help="Title for the HTML report")
def report(tree: Path, metadata: Optional[Path], logs: Optional[Path], title: str, outdir: Optional[Path], prefix: str, **context: Any):
    """
    Create phylogenetic tree from alignment

    Example usage:

        snippy-ng utils report --tree tree.newick
    """
    from snippy_ng import Context
    from snippy_ng.pipelines.report import ReportPipelineBuilder
    import json
    import csv

    if metadata:
        if metadata.suffix.lower() == ".json":
        # if metadata is a JSON file, read it and pass the content as a string to the pipeline
            with open(metadata, "r") as f:
                metadata = json.dumps(json.load(f))
        elif metadata.suffix.lower() == ".csv" or metadata.suffix.lower() == ".tsv":
            # if metadata is a CSV/TSV file, read it and convert to a list of dicts, then pass as a JSON string to the pipeline
            with open(metadata, "r") as f:
                if metadata.suffix.lower() == ".csv":
                    reader = csv.DictReader(f)
                else:
                    reader = csv.DictReader(f, delimiter="\t")
                metadata = json.dumps(list(reader))
        else:
            raise click.BadParameter("Metadata file must be in JSON, CSV, or TSV format", param_hint="--metadata")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = ReportPipelineBuilder(
        tree=tree,
        title=title,
        metadata=metadata,
        logs=logs,
        prefix=prefix,
    ).build()

    context["outdir"] = outdir or Path(".").absolute()
    # Run the pipeline
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)