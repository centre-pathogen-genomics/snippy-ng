import click
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, GlobalOption, add_snippy_global_options, create_outdir_callback
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=Path("tree"), required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=AbsolutePath), help="Output directory for phylogenetic tree results", callback=create_outdir_callback, cls=GlobalOption)
@click.option("--prefix", "-p", default="tree", help="Prefix for output files", cls=GlobalOption)
@add_snippy_global_options(exclude=['prefix', 'outdir'])
@click.argument("alignment", required=True, type=click.Path(exists=True, readable=True, path_type=AbsolutePath))
@click.option("--model", type=click.STRING, default="GTR+G4", help="Substitution model to use for tree construction. Use 'TEST' to select the best model.")
@click.option("--bootstrap", type=click.INT, default=1000, help="Number of bootstrap replicates to perform")
@click.option("--fconst", type=click.STRING, default=None, help="Constant sites counts (string or path to file)")
def tree(**config):
    """
    Create phylogenetic tree from alignment

    Example usage:

        snippy-ng tree core.full.aln
    """
    from snippy_ng.pipelines.tree import TreePipelineBuilder
    from snippy_ng.logging import logger

    #if fconst is a path read the content
    fconst = config.get("fconst")
    if not fconst and Path(config["alignment"]).with_suffix(".fconst").exists():
        logger.debug("No fconst provided, but found .fconst file corresponding to the alignment. Using this file for constant sites counts.")
        fconst = Path(config["alignment"]).with_suffix(".fconst")

    if fconst and Path(fconst).is_file():
        with open(fconst, 'r') as f:
            fconst = f.read().strip()
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = TreePipelineBuilder(
        aln=config["alignment"],
        model=config["model"],
        bootstrap=config["bootstrap"],
        fconst=fconst,
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
        prefix=config["prefix"],
    ).build()

    # Run the pipeline
    pipeline.run(
        skip_check=config['skip_check'],
        check=config['check'],
        cwd=config['outdir'],
        quiet=config['quiet'],
        create_missing=config['create_missing'],
        keep_incomplete=config['keep_incomplete'],
    )