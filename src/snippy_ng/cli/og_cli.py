from pathlib import Path
import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True}, short_help="Backwards-compatible SNP calling pipeline")
@snippy_global_options
@click.option("--reference", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--R1", "--pe1", "--left", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R2 (right)")
@click.option("--se", "--single", default='', type=click.STRING, help="Single-end reads")
@click.option("--ctgs", "--contigs", default='', type=click.STRING, help="Use these contigs instead of reads")
@click.option("--peil", default='', type=click.STRING, help="Paired-end interleaved reads")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--targets", default='', type=click.STRING, help="Only call SNPs from this BED file")
@click.option("--subsample", default=1.0, type=float, help="Subsample FASTQ to this proportion")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
@click.option("--report/--no-report", default=False, help="Produce report with visual alignment per variant")
@click.option("--cleanup/--no-cleanup", default=False, help="Remove unnecessary files (e.g., BAMs)")
@click.option("--rgid", default='', type=click.STRING, help="Use this @RG ID in the BAM header")
@click.option("--unmapped/--no-unmapped", default=False, help="Keep unmapped reads in BAM and write FASTQ")
@click.option("--mapqual", default=60, type=int, help="Minimum read mapping quality to consider")
@click.option("--basequal", default=13, type=int, help="Minimum base quality to consider")
@click.option("--mincov", default=10, type=int, help="Minimum site depth for calling alleles")
@click.option("--minfrac", default=0.0, type=float, help="Minimum proportion for variant evidence (0=AUTO)")
@click.option("--minqual", default=100.0, type=float, help="Minimum quality in VCF column 6")
@click.option("--maxsoft", default=10, type=int, help="Maximum soft clipping to allow")
@click.option("--bwaopt", default='', type=click.STRING, help="Extra BWA MEM options")
@click.option("--fbopt", default='', type=click.STRING, help="Extra Freebayes options")
def og(**kwargs):
    """
    Drop-in replacement for Snippy with feature parity.

    Examples:

        $ snippy-ng og --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.pipeline import Pipeline
    from snippy_ng.stages.setup import PrepareReference
    from snippy_ng.stages.alignment import BWAMEMReadsAligner, PreAlignedReads
    from snippy_ng.stages.calling import FreebayesCaller
    from snippy_ng.exceptions import DependencyError, MissingOutputError
    from snippy_ng.seq_utils import guess_format

    from pydantic import ValidationError

    def error(msg):
        click.echo(f"Error: {msg}", err=True)
        raise click.Abort()

    if not kwargs["force"] and kwargs["outdir"].exists():
        error(f"Output folder '{kwargs['outdir']}' already exists! Use --force to overwrite.")

    # check if output folder exists
    if not kwargs["outdir"].exists():
        kwargs["outdir"].mkdir(parents=True, exist_ok=True)

    # combine R1 and R2 into reads
    kwargs["reads"] = []
    if kwargs["r1"]:
        kwargs["reads"].append(kwargs["r1"])
    if kwargs["r2"]:
        kwargs["reads"].append(kwargs["r2"])
    if not kwargs["reads"] and not kwargs["bam"]:
        error("Please provide reads or a BAM file!")
    
    
    # Choose stages to include in the pipeline
    stages = []
    try:
        if Path(kwargs["reference"]).is_dir():
            # TODO use json file to get reference
            kwargs["reference"] = (Path(kwargs["reference"]) / "reference" / "ref.fa").resolve()
        else:
            reference_format = guess_format(kwargs["reference"])
            if not reference_format:
                error(f"Could not determine format of reference file '{kwargs['reference']}'")

            setup = PrepareReference(
                    input=kwargs["reference"],
                    ref_fmt=reference_format,
                    **kwargs,
                )
            kwargs["reference"] = setup.output.reference
            stages.append(setup)
        if kwargs["bam"]:
            aligner = PreAlignedReads(**kwargs)
        else:
            aligner = BWAMEMReadsAligner(**kwargs)
        kwargs["bam"] = aligner.output.bam
        stages.append(aligner)
        stages.append(FreebayesCaller(**kwargs))
    except ValidationError as e:
        error(e)
    
    # Move from CLI land into Pipeline land
    snippy = Pipeline(stages=stages)
    snippy.welcome()

    if not kwargs["skip_check"]:
        try:
            snippy.validate_dependencies()
        except DependencyError as e:
            snippy.error(f"Invalid dependencies! Please install '{e}' or use --skip-check to ignore.")
            return 1
    
    if kwargs["check"]:
        return 0

    # Set working directory to output folder
    snippy.set_working_directory(kwargs["outdir"])
    try:
        snippy.run(quiet=kwargs["quiet"])
    except MissingOutputError as e:
        snippy.error(e)
        return 1
    except RuntimeError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup()
    snippy.goodbye()

