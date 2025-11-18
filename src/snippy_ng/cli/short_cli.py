import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@snippy_global_options
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--R1", "--pe1", "--left", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R2 (right)")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--clean-reads", is_flag=True, default=False, help="Clean and filter reads with fastp before alignment")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2", "bwamem"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--freebayes-opts", default='', type=click.STRING, help="Extra options for Freebayes")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--min-depth", default=10, type=click.INT, help="Minimum coverage to call a variant")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
@click.option("--header", default=None, type=click.STRING, help="Header for the output FASTA file (if not provided, reference headers are kept)")
def short(**config):
    """
    Short read based SNP calling pipeline

    Examples:

        $ snippy-ng short --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.pipelines.short import create_short_stages
    import click
    from snippy_ng.snippy import Snippy
    from snippy_ng.exceptions import DependencyError, MissingOutputError
    
    # combine R1 and R2 into reads
    config["reads"] = []
    
    if config.get("r1"):
        config["reads"].append(config["r1"])
    if config.get("r2"):
        config["reads"].append(config["r2"])
    if not config["reads"] and not config.get("bam"):
        click.UsageError("Please provide reads or a BAM file!")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    stages = create_short_stages(config)
    
    # Run the pipeline
    # Move from CLI land into Pipeline land
    snippy = Snippy(stages=stages)
    snippy.welcome()

    if not config.get("skip_check", False):
        try:
            snippy.validate_dependencies()
        except DependencyError as e:
            snippy.error(f"Invalid dependencies! Please install '{e}' or use --skip-check to ignore.")
            raise click.Abort()
    
    if config.get("check", False):
        return 0

    # Set working directory to output folder
    snippy.set_working_directory(config["outdir"])
    try:
        snippy.run(
            quiet=config.get("quiet", False),
            continue_last_run=config.get("continue", False),
            keep_incomplete=config.get("keep_incomplete", False),
        )
    except MissingOutputError as e:
        snippy.error(e)
        raise click.Abort()
    except RuntimeError as e:
        snippy.error(e)
        raise click.Abort()
    
    snippy.cleanup()
    snippy.goodbye()
    
