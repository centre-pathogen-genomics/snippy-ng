import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank) or prepared reference directory")
@click.option("--reads", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Long reads file (FASTQ)")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--clean-reads", is_flag=True, default=True, help="Remove short and low-quality reads before alignment")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--minimap-preset", default="map-ont", type=click.Choice(["map-ont", "lr:hq", "map-hifi", "map-pb"]), help="Preset for minimap2 alignment")
@click.option("--caller", default="clair3", type=click.Choice(["clair3", "freebayes"]), help="Variant caller to use")
@click.option("--caller-opts", default="", type=click.STRING, help="Additional options to pass to the variant caller")
@click.option("--clair3-model", default=None, type=click.Path(resolve_path=True), help="Absolute path to Clair3 model file.")
@click.option("--clair3-fast-mode", is_flag=True, default=False, help="Enable fast mode in Clair3 for quicker variant calling")
@click.option("--min-read-len", type=click.INT, default=1000, help="Minimum read length to keep when cleaning reads")
@click.option("--min-read-qual", type=click.FLOAT, default=10, help="Minimum read quality to keep when cleaning reads")
@click.option("--min-depth", default=10, type=click.INT, help="Minimum coverage to call a variant")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
def long(**config):
    """
    Long read based SNP calling pipeline

    Examples:

        $ snippy-ng long --reference ref.fa --reads long_reads.fq --outdir output
    """
    from snippy_ng.pipelines.long import LongPipelineBuilder
    import click
    
    if not config.get("reads") and not config.get("bam"):
        raise click.UsageError("Please provide reads or a BAM file!")

    if config.get("caller") == "clair3" and not config.get("clair3_model"):
        raise click.UsageError("Please provide a Clair3 model file (--clair3-model) when using Clair3 as the variant caller!")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = LongPipelineBuilder(
        reference=config["reference"],
        reads=config["reads"],
        prefix=config["prefix"],
        bam=config["bam"],
        aligner=config["aligner"],
        aligner_opts=config["aligner_opts"],
        minimap_preset=config["minimap_preset"],
        caller=config["caller"],
        caller_opts=config["caller_opts"],
        clair3_model=config.get("clair3_model"),
        clair3_fast_mode=config["clair3_fast_mode"],
        downsample=config["downsample"],
        min_read_len=config.get("min_read_len"),
        min_read_qual=config.get("min_read_qual"),
        min_qual=config["min_qual"],
        min_depth=config["min_depth"],
        mask=config["mask"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    ).build()
    
    # Run the pipeline
    return pipeline.run(
        skip_check=config['skip_check'],
        check=config['check'],
        cwd=config['outdir'],
        quiet=config['quiet'],
        create_missing=config['create_missing'],
        keep_incomplete=config['keep_incomplete'],
    )
    
