import click
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options, create_outdir_callback, GlobalOption
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=Path("reference"), required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=AbsolutePath), help="Output directory for the prepared reference", callback=create_outdir_callback, cls=GlobalOption)
@add_snippy_global_options(exclude=['outdir', 'prefix'])
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, readable=True, path_type=AbsolutePath), help="Reference genome (FASTA or GenBank)")
def ref(**config):
    """
    Utility to prepare a reference genome for use with snippy-ng. 
    
    This includes indexing the reference and creating any necessary auxiliary files.

    Examples:

        $ snippy-ng ref --reference ref.fa --outdir output
    """
    from snippy_ng.pipelines.common import prepare_reference
    from snippy_ng.pipelines import SnippyPipeline

    ref_stage = prepare_reference(
                reference_path=config["reference"],
                output_directory=config["outdir"]
            )
    pipeline = SnippyPipeline(stages=[ref_stage])

    pipeline.run(
        skip_check=config["skip_check"],
        check=config["check"],
        cwd=Path("."),
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )