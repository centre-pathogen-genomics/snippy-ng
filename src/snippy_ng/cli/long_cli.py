import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@snippy_global_options
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--reads", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Long reads file (FASTQ)")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--clair3-model", default=None, type=click.Path(resolve_path=True), help="Path to Clair3 model file. If not provided, freebayes will be used for variant calling.")
@click.option("--clair3-fast-mode", is_flag=True, default=False, help="Enable fast mode in Clair3 for quicker variant calling")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--clean-reads", is_flag=True, default=True, help="Remove short and low-quality reads before alignment")
@click.option("--min-read-len", type=click.INT, default=1000, help="Minimum read length to keep when cleaning reads")
@click.option("--min-read-qual", type=click.FLOAT, default=10, help="Minimum read quality to keep when cleaning reads")
@click.option("--min-depth", default=10, type=click.INT, help="Minimum coverage to call a variant")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
def long(**config):
    """
    Long read based SNP calling pipeline

    Examples:

        $ snippy-ng long --reference ref.fa --reads long_reads.fq --outdir output
    """
    from snippy_ng.pipelines.long import create_long_pipeline_stages
    from snippy_ng.cli.utils.pipeline_runner import run_snippy_pipeline
    import click
    
    if not config.get("reads") and not config.get("bam"):
        raise click.UsageError("Please provide reads or a BAM file!")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    stages = create_long_pipeline_stages(
        reference=config["reference"],
        reads=config["reads"],
        prefix=config["prefix"],
        bam=config["bam"],
        clair3_model=config.get("clair3_model"),
        downsample=config["downsample"],
        min_read_len=config.get("min_read_len"),
        min_read_qual=config.get("min_read_qual"),
        min_qual=config["min_qual"],
        min_depth=config["min_depth"],
        mask=config["mask"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    )
    
    # Run the pipeline
    return run_snippy_pipeline(config, stages)
    
