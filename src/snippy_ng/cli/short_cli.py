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
def short(**config):
    """
    Short read based SNP calling pipeline

    Examples:

        $ snippy-ng short --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.pipelines.short import create_short_pipeline_stages
    from snippy_ng.cli.utils.pipeline_runner import run_snippy_pipeline
    import click
    
    # combine R1 and R2 into reads
    reads = []
    if config.get("r1"):
        reads.append(config["r1"])
    if config.get("r2"):
        reads.append(config["r2"])
    if not reads and not config.get("bam"):
        raise click.UsageError("Please provide reads or a BAM file!")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    stages = create_short_pipeline_stages(
        reference=config["reference"],
        reads=reads,
        prefix=config["prefix"],
        bam=config["bam"],
        clean_reads=config["clean_reads"],
        downsample=config["downsample"],
        aligner=config["aligner"],
        aligner_opts=config["aligner_opts"],
        freebayes_opts=config["freebayes_opts"],
        mask=config["mask"],
        min_depth=config["min_depth"],
        min_qual=config["min_qual"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    )
    
    # Run the pipeline
    return run_snippy_pipeline(config, stages)
    
