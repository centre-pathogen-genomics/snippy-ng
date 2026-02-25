import click
from typing import Any, Optional
from pathlib import Path
from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=AbsolutePath(exists=True, readable=True), help="Reference genome (FASTA or GenBank) or prepared reference directory")
@click.option("--R1", "--pe1", "--left", default=None, type=AbsolutePath(exists=True, readable=True), help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default=None, type=AbsolutePath(exists=True, readable=True), help="Reads, paired-end R2 (right)")
@click.option("--bam", default=None, type=AbsolutePath(exists=True), help="Use this BAM file instead of aligning reads")
@click.option("--clean-reads", is_flag=True, default=False, help="Clean and filter reads with fastp before alignment")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--mask", default=None, type=AbsolutePath(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2", "bwamem"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--caller-opts", default='', type=click.STRING, help="Extra options for Freebayes")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
def short(reference: Path, r1: Optional[Path], r2: Optional[Path], bam: Optional[Path], clean_reads: bool, downsample: Optional[float], mask: Optional[Path], aligner: str, aligner_opts: str, caller_opts: str, min_qual: float, prefix: str, outdir: Path, **context: Any):
    """
    Short read based SNP calling pipeline

    Examples:

        $ snippy-ng short --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.context import Context
    from snippy_ng.pipelines.short import ShortPipelineBuilder
    
    # combine R1 and R2 into reads
    reads = []
    if r1:
        reads.append(r1)
    if r2:
        reads.append(r2)
    if not reads and not bam:
        raise click.UsageError("Please provide reads or a BAM file!")
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = ShortPipelineBuilder(
        reference=reference,
        reads=reads,
        prefix=prefix,
        bam=bam,
        clean_reads=clean_reads,
        downsample=downsample,
        aligner=aligner,
        aligner_opts=aligner_opts,
        caller_opts=caller_opts,
        mask=mask,
        min_qual=min_qual,
    ).build()
    
    # Run the pipeline
    context["outdir"] = outdir
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
    
