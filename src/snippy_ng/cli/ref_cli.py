import click
from typing import Any
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options, check_outdir_callback, GlobalOption
from pathlib import Path


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@click.option("--outdir", "-o", default=Path("reference"), required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True), help="Output directory for the prepared reference", callback=check_outdir_callback, cls=GlobalOption)
@add_snippy_global_options(exclude=['outdir', 'prefix'])
@click.argument("reference", required=True, type=AbsolutePath(exists=True, readable=True))
def ref(reference: Path, **context: Any):
    """
    Prepare a reference genome for use with snippy-ng
    
    This includes converting gbk to gff, indexing the reference and creating any necessary auxiliary files.

    Examples:

        $ snippy-ng utils ref my_reference.gbk
    """
    from snippy_ng.context import Context
    from snippy_ng.pipelines.common import prepare_reference
    from snippy_ng.pipelines import SnippyPipeline

    ref_stage = prepare_reference(
                reference_path=reference,
            )
    pipeline = SnippyPipeline(stages=[ref_stage])

    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
