import click
from typing import Any, Optional
from pathlib import Path
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=AbsolutePath(exists=True, readable=True), help="Reference genome (FASTA or GenBank) or prepared reference directory")
@click.option("--reads", default=None, type=AbsolutePath(exists=True, readable=True), help="Long reads file (FASTQ)")
@click.option("--bam", default=None, type=AbsolutePath(exists=True), help="Use this BAM file instead of aligning reads")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--clean-reads", is_flag=True, default=True, help="Remove short and low-quality reads before alignment")
@click.option("--mask", default=None, type=AbsolutePath(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--minimap-preset", default="map-ont", type=click.Choice(["map-ont", "lr:hq", "map-hifi", "map-pb"]), help="Preset for minimap2 alignment")
@click.option("--caller", default="clair3", type=click.Choice(["clair3", "freebayes"]), help="Variant caller to use")
@click.option("--caller-opts", default="", type=click.STRING, help="Additional options to pass to the variant caller")
@click.option("--clair3-model", default=None, type=AbsolutePath(), help="Absolute path to Clair3 model file.")
@click.option("--clair3-fast-mode", is_flag=True, default=False, help="Enable fast mode in Clair3 for quicker variant calling")
@click.option("--min-read-len", type=click.INT, default=1000, help="Minimum read length to keep when cleaning reads")
@click.option("--min-read-qual", type=click.FLOAT, default=10, help="Minimum read quality to keep when cleaning reads")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
def long(reference: Path, reads: Optional[Path], bam: Optional[Path], downsample: Optional[float], clean_reads: bool, mask: Optional[Path], aligner: str, aligner_opts: str, minimap_preset: str, caller: str, caller_opts: str, clair3_model: Optional[Path], clair3_fast_mode: bool, min_read_len: int, min_read_qual: float, min_qual: float, prefix: str, outdir: Path, **context: Any):
    """
    Long read based SNP calling pipeline

    Examples:

        $ snippy-ng long --reference ref.fa --reads long_reads.fq --outdir output
    """
    from snippy_ng import Context
    from snippy_ng.pipelines.long import LongPipelineBuilder
    import click
    
    if not reads and not bam:
        raise click.UsageError("Please provide reads or a BAM file!")

    if caller == "clair3" and not clair3_model:
        raise click.UsageError("Please provide a Clair3 model file (--clair3-model) when using Clair3 as the variant caller!")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = LongPipelineBuilder(
        reference=reference,
        reads=reads,
        prefix=prefix,
        bam=bam,
        aligner=aligner,
        aligner_opts=aligner_opts,
        minimap_preset=minimap_preset,
        caller=caller,
        caller_opts=caller_opts,
        clair3_model=clair3_model,
        clair3_fast_mode=clair3_fast_mode,
        downsample=downsample,
        min_read_len=min_read_len,
        min_read_qual=min_read_qual,
        min_qual=min_qual,
        mask=mask,
    ).build()
    
    # Run the pipeline
    context["outdir"] = outdir
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
    
