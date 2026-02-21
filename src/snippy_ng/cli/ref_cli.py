import click
from typing import Any
from snippy_ng.cli.utils import absolute_path_callback
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options, create_outdir_callback, GlobalOption
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=Path("reference"), required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True), help="Output directory for the prepared reference", callback=create_outdir_callback, cls=GlobalOption)
@add_snippy_global_options(exclude=['outdir', 'prefix'])
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Reference genome (FASTA or GenBank)")
def ref(reference: Path, outdir: Path, **context: Any):
    """
    Prepare a reference genome for use with snippy-ng. 
    
    This includes converting gbk to gff, indexing the reference and creating any necessary auxiliary files.

    Examples:

        $ snippy-ng utils ref --reference ref.fa --outdir output
    """
    from snippy_ng import Context
    from snippy_ng.pipelines.common import prepare_reference
    from snippy_ng.pipelines import SnippyPipeline

    ref_stage = prepare_reference(
                reference_path=reference,
                output_directory=outdir
            )
    pipeline = SnippyPipeline(stages=[ref_stage])

    context["outdir"] = Path(".")
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)