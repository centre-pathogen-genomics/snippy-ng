import click
from typing import Any, Optional, Literal
from pathlib import Path
from snippy_ng.cli.utils import AbsolutePath, reference_or_accession_callback
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=click.STRING, callback=reference_or_accession_callback, help="Reference genome (FASTA or GenBank), prepared reference directory, or NCBI GCF/GCA assembly accession")
@click.option("--reads", default=None, type=AbsolutePath(exists=True, readable=True), help="Long reads file (FASTQ)")
@click.option("--bam", default=None, type=AbsolutePath(exists=True), help="Use this BAM file instead of aligning reads")
@click.option("--clean-reads/--no-clean-reads", is_flag=True, default=True, help="Remove short and low-quality reads before alignment")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--min-read-len", type=click.INT, default=1000, help="Minimum read length to keep when cleaning reads")
@click.option("--min-read-qual", type=click.FLOAT, default=10, help="Minimum read quality to keep when cleaning reads")
@click.option("--mask", default=None, type=AbsolutePath(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--depth-mask", default=0, type=click.INT, help="Mask regions in the output fasta with Ns if the read depth is below this threshold")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--minimap-preset", default="lr:hq", type=click.Choice(["lr:hq", "map-ont", "map-hifi", "map-pb"]), help="Preset for minimap2 alignment")
@click.option("--caller", default="clair3", type=click.Choice(["clair3", "freebayes"]), help="Variant caller to use")
@click.option("--caller-opts", default="", type=click.STRING, help="Additional options to pass to the variant caller")
@click.option("--caller-map-qual", default=60, type=click.INT, help="Minimum mapping quality for caller to consider a read")
@click.option("--clair3-model", default=None, type=AbsolutePath(), help="Path to Clair3 model file. If not provided, will attempt to find a suitable model using LongBow")
@click.option("--min-qual", default=None, type=click.FLOAT, help="Minimum QUAL threshold for low quality variant masking. Default is AUTO for Clair3 and 100 for FreeBayes")
@click.option("--report/--no-report", default=False, help="Create a per-sample HTML report")
def long(
    reference: Path | str,
    reads: Optional[Path],
    bam: Optional[Path],
    downsample: Optional[float],
    clean_reads: bool,
    min_read_len: int,
    min_read_qual: float,
    mask: Optional[Path],
    depth_mask: int,
    aligner: str,
    aligner_opts: str,
    minimap_preset: str,
    caller: Literal["clair3", "freebayes"],
    caller_opts: str,
    caller_map_qual: int,
    clair3_model: Optional[Path],
    min_qual: Optional[float],
    report: bool,
    prefix: str,
    **context: Any,
):
    """
    Long read based SNP calling pipeline

    Examples:

        $ snippy-ng long --reference ref.fa --reads long_reads.fq --outdir output
    """
    from snippy_ng.context import Context
    from snippy_ng.pipelines.long import LongPipelineBuilder
    import click
    
    if not reads and not bam:
        raise click.UsageError("Please provide reads or a BAM file!")

    if caller == "clair3" and not clair3_model and not reads:
        raise click.UsageError("Please provide --clair3-model when using Clair3 with BAM/CRAM input only.")
    
    # Convert reference to accession if it's a string, otherwise keep as Path
    reference_accession = reference if isinstance(reference, str) else None
    reference_path = None if reference_accession else Path(reference)
    
    if min_qual is None and caller == "freebayes":
        min_qual = 100.0
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = LongPipelineBuilder(
        reference=reference_path,
        reference_accession=reference_accession,
        reads=reads,
        prefix=prefix,
        bam=bam,
        aligner=aligner,
        aligner_opts=aligner_opts,
        minimap_preset=minimap_preset,
        caller=caller,
        caller_opts=caller_opts,
        clair3_model=clair3_model,
        clean_reads=clean_reads,
        downsample=downsample,
        min_read_len=min_read_len,
        min_read_qual=min_read_qual,
        min_qual=min_qual,
        min_mapping_quality=caller_map_qual,
        mask=mask,
        depth_mask=depth_mask,
        report=report,
    ).build()
    
    # Run the pipeline
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
    
