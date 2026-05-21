import click
from typing import Any, Optional
from pathlib import Path
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options

@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=AbsolutePath(exists=True, resolve_path=False, readable=True), help="Reference genome (FASTA or GenBank) or prepared reference directory")
@click.option("--assembly", "--asm", required=True, type=AbsolutePath(exists=True, readable=True), help="Assembly in FASTA format")
@click.option("--mask", default=None, type=AbsolutePath(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--caller", default="nucmer", type=click.Choice(["nucmer", "paftools"]), help="Caller to use for assembly-based SNP calling")
@click.option("--caller-opts", default="", type=click.STRING, help="Extra options for the assembly caller")
@click.option("--report/--no-report", default=False, help="Create a per-sample HTML report")
def asm(reference: Path, assembly: Path, mask: Optional[Path], caller: str, caller_opts: str, report: bool, prefix: str, **context: Any):
    """
    Assembly based SNP calling pipeline

    Examples:

s    """
    from snippy_ng.context import Context
    from snippy_ng.pipelines.asm import AsmPipelineBuilder
    
    # Choose stages to include in the pipeline
    # ensure this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = AsmPipelineBuilder(
        reference=reference,
        assembly=assembly,
        prefix=prefix,
        caller=caller,
        caller_opts=caller_opts,
        mask=mask,
        report=report,
    ).build()

    # Create a context object to pass to the pipeline run method
    run_ctx = Context(**context)
    # Run the pipeline with the provided context
    return pipeline.run(run_ctx)
