import click
from typing import Any, Optional
from pathlib import Path
from snippy_ng.cli.utils import absolute_path_callback
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options

@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=False, readable=True), callback=absolute_path_callback, help="Reference genome (FASTA or GenBank) or prepared reference directory")
@click.option("--assembly", "--asm", required=True, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Assembly in FASTA format")
@click.option("--mask", default=None, type=click.Path(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
def asm(reference: Path, assembly: Path, mask: Optional[Path], prefix: str, **context: Any):
    """
    Assembly based SNP calling pipeline

    Examples:

        $ snippy-ng asm --reference ref.fa --assembly assembly.fa --outdir output
    """
    from snippy_ng import Context
    from snippy_ng.pipelines.asm import AsmPipelineBuilder
    
    # Choose stages to include in the pipeline
    # ensure this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = AsmPipelineBuilder(
        reference=reference,
        assembly=assembly,
        prefix=prefix,
        mask=mask,
    ).build()

    # Create a context object to pass to the pipeline run method
    run_ctx = Context(**context)
    # Run the pipeline with the provided context
    return pipeline.run(run_ctx)

