import click
from snippy_ng.cli.utils import absolute_path_callback
from snippy_ng.cli.utils.globals import CommandWithGlobals, GlobalOption, add_snippy_global_options, create_outdir_callback
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=None, required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True), help="Output directory for phylogenetic tree results", callback=create_outdir_callback, cls=GlobalOption)
@click.option("--prefix", "-p", default="report", help="Prefix for html report", cls=GlobalOption)
@add_snippy_global_options(exclude=['prefix', 'outdir'])
@click.option("--tree", required=True, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Newick tree file to include in the report")
@click.option("--metadata", required=False, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Optional metadata file (JSON or CSV) to include in the report")
@click.option("--logs", required=False, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Optional log file to include in the report")
@click.option("--title", required=False, type=click.STRING, default="Snippy-NG Report", help="Title for the HTML report")
def report(**config):
    """
    Create phylogenetic tree from alignment

    Example usage:

        snippy-ng utils report --tree tree.newick
    """
    from snippy_ng.pipelines.report import ReportPipelineBuilder
    import json
    import csv

    if config.get("metadata"):
        if config["metadata"].suffix.lower() == ".json":
        # if metadata is a JSON file, read it and pass the content as a string to the pipeline
            with open(config["metadata"], "r") as f:
                config["metadata"] = json.dumps(json.load(f))
        elif config["metadata"].suffix.lower() == ".csv" or config["metadata"].suffix.lower() == ".tsv":
            # if metadata is a CSV/TSV file, read it and convert to a list of dicts, then pass as a JSON string to the pipeline
            with open(config["metadata"], "r") as f:
                if config["metadata"].suffix.lower() == ".csv":
                    reader = csv.DictReader(f)
                else:
                    reader = csv.DictReader(f, delimiter="\t")
                config["metadata"] = json.dumps(list(reader))
        else:
            raise click.BadParameter("Metadata file must be in JSON, CSV, or TSV format", param_hint="--metadata")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = ReportPipelineBuilder(
        tree=config["tree"],
        title=config["title"],
        metadata=config.get("metadata"),
        logs=config.get("logs"),
        prefix=config["prefix"],
    ).build()

    outdir = config["outdir"] or Path(".").absolute()
    # Run the pipeline
    pipeline.run(
        skip_check=config['skip_check'],
        check=config['check'],
        outdir=outdir,
        quiet=config['quiet'],
        create_missing=config['create_missing'],
        keep_incomplete=config['keep_incomplete'],
    )