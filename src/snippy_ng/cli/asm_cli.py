import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@snippy_global_options
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--assembly", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Assembly in FASTA format")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
def asm(**config):
    """
    Assembly based SNP calling pipeline

    Examples:

        $ snippy-ng asm --reference ref.fa --assembly assembly.fa --outdir output
    """
    from snippy_ng.pipelines.asm import create_asm_pipeline_stages
    from snippy_ng.cli.utils.pipeline_runner import run_snippy_pipeline
    
    # Choose stages to include in the pipeline
    # ensure this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    stages = create_asm_pipeline_stages(
        reference=config["reference"],
        assembly=config["assembly"],
        prefix=config["prefix"],
        mask=config["mask"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    )
    
    # Run the pipeline
    return run_snippy_pipeline(config, stages)
