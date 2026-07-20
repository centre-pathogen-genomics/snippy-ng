import click
from typing import Any, Optional, Literal
from pathlib import Path
from snippy_ng.cli.utils import AbsolutePath, assembly_or_accession_callback, reference_or_accession_callback, resolve_cli_input
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options

@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=click.STRING, callback=reference_or_accession_callback, help="Reference genome (FASTA or GenBank), prepared reference directory, or NCBI GCF/GCA assembly accession")
@click.option("--assembly", "--asm", default=None, type=click.STRING, callback=assembly_or_accession_callback, help="Assembly in FASTA format or NCBI GCF/GCA or AllTheBacteria SAMN/SAMEA accession")
@click.argument("assembly_arg", required=False, metavar="ASSEMBLY", type=click.STRING, callback=assembly_or_accession_callback)
@click.option("--vcf", default=None, type=AbsolutePath(exists=True), help="Use this VCF file instead of calling variants")
@click.option("--mask", default=None, type=AbsolutePath(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--caller", default="nucmer", type=click.Choice(["nucmer", "paftools"]), help="Caller to use for assembly-based SNP calling")
@click.option("--caller-opts", default="", type=click.STRING, help="Extra options for the assembly caller")
@click.option("--report/--no-report", default=False, help="Create a per-sample HTML report")
def asm(reference: Path | str, assembly: Optional[Path | str], assembly_arg: Optional[Path | str], vcf: Optional[Path], mask: Optional[Path], caller: Literal["nucmer", "paftools"], caller_opts: str, report: bool, prefix: str, sample_name: Optional[str], **context: Any):
    """
    Assembly based SNP calling pipeline

    Examples:
    
        $ snippy-ng asm --reference ref.fa asm.fa --outdir output
    
    """
    from snippy_ng.context import Context
    from snippy_ng.pipelines.asm import AsmPipelineBuilder

    assembly = resolve_cli_input(
        assembly,
        assembly_arg,
        option_name="--assembly/--asm",
        arg_name="assembly",
    )
    if assembly is None:
        raise click.UsageError("Please provide an assembly via the positional argument or --assembly/--asm.")

    # Convert reference to accession if it's a string, otherwise keep as Path
    reference_accession = reference if isinstance(reference, str) else None
    reference_path = None if reference_accession else Path(reference)
    assembly_accession = assembly if isinstance(assembly, str) else None
    assembly_path = None if assembly_accession else Path(assembly)
    
    # Choose stages to include in the pipeline
    # ensure this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = AsmPipelineBuilder(
        reference=reference_path,
        reference_accession=reference_accession,
        assembly=assembly_path,
        assembly_accession=assembly_accession,
        vcf=vcf,
        prefix=prefix,
        sample_name=sample_name,
        caller=caller,
        caller_opts=caller_opts,
        mask=mask,
        report=report,
    ).build()

    # Create a context object to pass to the pipeline run method
    run_ctx = Context(**context)
    # Run the pipeline with the provided context
    return pipeline.run(run_ctx)
