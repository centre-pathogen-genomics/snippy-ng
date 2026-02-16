import click
from snippy_ng.cli.utils import absolute_path_callback
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options

@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=False, readable=True), callback=absolute_path_callback, help="Reference genome (FASTA or GenBank) or prepared reference directory")
@click.option("--assembly", "--asm", required=True, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Assembly in FASTA format")
@click.option("--mask", default=None, type=click.Path(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
def asm(**config):
    """
    Assembly based SNP calling pipeline

    Examples:

        $ snippy-ng asm --reference ref.fa --assembly assembly.fa --outdir output
    """
    from snippy_ng.pipelines.asm import AsmPipelineBuilder
    
    # Choose stages to include in the pipeline
    # ensure this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = AsmPipelineBuilder(
        reference=config["reference"],
        assembly=config["assembly"],
        prefix=config["prefix"],
        mask=config["mask"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
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

