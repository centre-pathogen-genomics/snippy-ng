import click
from typing import Any, Optional
from pathlib import Path
from snippy_ng.cli.utils import AbsolutePath, reference_or_accession_callback, reads_or_accession_callback, resolve_cli_input, is_sra_accession
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@add_snippy_global_options()
@click.option("--reference", "--ref", required=True, type=click.STRING, callback=reference_or_accession_callback, help="Reference genome (FASTA or GenBank), prepared reference directory, or NCBI GCF/GCA assembly accession")
@click.option("--R1", "--pe1", "--left", default=None, type=click.STRING, callback=reads_or_accession_callback, help="Reads, paired-end R1 (left) or SRA accession (SRR/ERR/DRR)")
@click.option("--R2", "--pe2", "--right", default=None, type=AbsolutePath(exists=True, readable=True), help="Reads, paired-end R2 (right)")
@click.argument("read_args", nargs=-1, metavar="READS", type=click.STRING, callback=reads_or_accession_callback)
@click.option("--bam", default=None, type=AbsolutePath(exists=True), help="Use this BAM file instead of aligning reads")
@click.option("--vcf", default=None, type=AbsolutePath(exists=True), help="Use this VCF file instead of calling variants")
@click.option("--clean-reads/--no-clean-reads", is_flag=True, default=False, help="Clean and filter reads with fastp before alignment")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--min-read-len", type=click.INT, default=15, help="Minimum read length to keep when cleaning reads")
@click.option("--min-read-qual", type=click.FLOAT, default=20, help="Minimum read quality to keep when cleaning reads")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2", "bwamem"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--caller-opts", default='', type=click.STRING, help="Extra options for Freebayes")
@click.option("--caller-map-qual", default=60, type=click.INT, help="Minimum mapping quality for caller to consider a read")
@click.option("--mask", default=None, type=AbsolutePath(exists=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--depth-mask", default=0, type=click.INT, help="Mask regions in the output fasta with Ns if the read depth is below this threshold")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Mark variants below this QUAL threshold as LowQual in the output VCF")
@click.option("--report/--no-report", default=False, help="Create a per-sample HTML report")
def short(
    reference: Path | str,
    r1: Optional[Path | str],
    r2: Optional[Path],
    read_args: tuple[Path | str, ...],
    bam: Optional[Path],
    vcf: Optional[Path],
    downsample: Optional[float],
    clean_reads: bool,
    min_read_len: int,
    min_read_qual: float,
    mask: Optional[Path],
    depth_mask: int,
    aligner: str,
    aligner_opts: str,
    caller_opts: str,
    caller_map_qual: int,
    min_qual: float,
    report: bool,
    prefix: str,
    sample_name: Optional[str],
    **context: Any,
):
    """
    Short read based SNP calling pipeline

    Examples:

        $ snippy-ng short --reference ref.fa reads_1.fq reads_2.fq --outdir output
        $ snippy-ng short --reference ref.fa SRR1234567 --outdir output
    """
    from snippy_ng.context import Context
    from snippy_ng.pipelines.short import ShortPipelineBuilder

    if len(read_args) > 2:
        raise click.UsageError("Please provide at most two positional read files: R1 [R2].")

    arg_r1 = read_args[0] if len(read_args) >= 1 else None
    arg_r2 = read_args[1] if len(read_args) >= 2 else None
    
    # Detect if arg_r1 or --R1 is a read accession
    read_accession = None

    # Check positional arg first
    if arg_r1 and isinstance(arg_r1, str) and is_sra_accession(arg_r1):
        read_accession = arg_r1
        arg_r1 = None
    # Check --R1 option
    elif r1 and isinstance(r1, str) and is_sra_accession(r1):
        read_accession = r1
        r1 = None

    r1 = resolve_cli_input(r1, arg_r1, option_name="--R1/--pe1/--left", arg_name="R1")
    r2 = resolve_cli_input(r2, arg_r2, option_name="--R2/--pe2/--right", arg_name="R2")

    reads = []
    if r1:
        reads.append(r1 if isinstance(r1, Path) else Path(r1))
    if r2:
        reads.append(r2)
    
    if not reads and not read_accession and not bam:
        raise click.UsageError("Please provide reads, a read accession, or a BAM file!")

    if vcf and not bam:
        raise click.UsageError("Please provide --bam when using --vcf; the alignment is required for depth masks and QC.")
    
    # Convert reference to accession if it's a string, otherwise keep as Path
    reference_accession = reference if isinstance(reference, str) else None
    reference_path = None if reference_accession else Path(reference)
    
    # Choose stages to include in the pipeline
    # this will raise ValidationError if config is invalid
    # we let this happen as we want to catch all config errors
    # before starting the pipeline
    pipeline = ShortPipelineBuilder(
        reference=reference_path,
        reference_accession=reference_accession,
        reads=reads,
        read_accession=read_accession,
        prefix=prefix,
        sample_name=sample_name,
        bam=bam,
        vcf=vcf,
        clean_reads=clean_reads,
        min_read_len=min_read_len,
        min_read_qual=min_read_qual,
        downsample=downsample,
        aligner=aligner,
        aligner_opts=aligner_opts,
        caller_opts=caller_opts,
        mask=mask,
        depth_mask=depth_mask,
        min_qual=min_qual,
        min_mapping_quality=caller_map_qual,
        report=report,
    ).build()
    
    # Run the pipeline
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
