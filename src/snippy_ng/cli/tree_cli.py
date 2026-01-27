import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options(exclude=['prefix', 'outdir'])
@click.argument("alignment", required=True, type=click.Path(exists=True, resolve_path=True, readable=True))
@click.option("--model", type=click.STRING, default="GTR+G4", help="Substitution model to use for tree construction")
@click.option("--bootstrap", type=click.INT, default=1000, help="Number of bootstrap replicates to perform")
@click.option("--fconst", type=click.STRING, default=None, help="Constant sites counts (string or path to file)")
def tree(**config):
    """
    Create phylogenetic tree from alignment

    Example usage:

        snippy-ng tree core.full.aln
    """
    from snippy_ng.pipelines.tree import create_tree_pipeline_stages
    from snippy_ng.cli.utils.pipeline_runner import run_snippy_pipeline

    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    stages = create_tree_pipeline_stages(
        aln=config["alignment"],
        model=config["model"],
        bootstrap=config["bootstrap"],
        fconst=config["fconst"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    )

    # Run the pipeline
    return run_snippy_pipeline(
        stages,
        skip_check=config['skip_check'],
        check=config['check'],
        outdir=Path('.'),
        quiet=config['quiet'],
        continue_last_run=config['continue_last_run'],
        keep_incomplete=config['keep_incomplete'],
    )