import click

from snippy_ng.__about__ import __version__, EXE, AUTHOR, URL
from pathlib import Path

def show_citation(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"Please cite {EXE} in your research: ...")
    ctx.exit()

def version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"{EXE} version {__version__}")
    ctx.exit()

@click.group(context_settings={"help_option_names": ["-h", "--help"]}, invoke_without_command=True)
@click.version_option(version=__version__, prog_name=EXE)
@click.option("--citation", is_flag=True, callback=show_citation, expose_value=False, help="Print citation for referencing this tool")
def snippy_ng():
    """
    Snippy-NG: The Next Generation of Variant Calling.
    """

@snippy_ng.command()
@click.option("--check/--no-check", default=False, help="Check dependencies are installed then exit")
@click.option("--skip-check/--no-skip-check", default=False, help="Skip dependency checks")
@click.option("--force/--no-force", default=False, help="Force overwrite of existing output folder")
@click.option("-q", "--quiet/--no-quiet", default=False, help="No screen output")
@click.option("--cpus", default=8, type=int, help="Maximum number of CPU cores to use")
@click.option("--ram", default=8, type=int, help="Try and keep RAM under this many GB")
@click.option("--tmpdir", default='/tmp', type=click.Path(), help="Fast temporary storage eg. local SSD")
@click.option("--reference", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA, GenBank, EMBL)")
@click.option("--R1", "--pe1", "--left", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R2 (right)")
@click.option("--se", "--single", default='', type=click.STRING, help="Single-end reads")
@click.option("--ctgs", "--contigs", default='', type=click.STRING, help="Use these contigs instead of reads")
@click.option("--peil", default='', type=click.STRING, help="Paired-end interleaved reads")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--targets", default='', type=click.STRING, help="Only call SNPs from this BED file")
@click.option("--subsample", default=1.0, type=float, help="Subsample FASTQ to this proportion")
@click.option("--outdir", required=True, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=Path), help="Output folder")
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
def snps(**kwargs):
    """
    Drop-in replacement for Snippy with feature parity.

    Examples:

        $ snippy-ng snps --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
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
    
    reference_format = guess_format(kwargs["reference"])
    if not reference_format:
        error(f"Could not determine format of reference file '{kwargs['reference']}'")

    # Choose stages to include in the pipeline
    stages = []
    try:
        setup = PrepareReference(input=kwargs["reference"], ref_fmt=reference_format)
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
        error_msg = "\n".join([f"{e['loc'][0]} ({e['input']}). {e['msg']}." for e in e.errors()])
        error(error_msg)
    
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




