from pathlib import Path
import click
from snippy_ng.cli.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True}, short_help="Run SNP calling pipeline")
@snippy_global_options
@click.option("--reference", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA, GenBank, EMBL)")
@click.option("--R1", "--pe1", "--left", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R2 (right)")
@click.option("--clean-reads", is_flag=True, default=False, help="Clean and filter reads with fastp before alignment")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2", "bwamem"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
def short(**kwargs):
    """
    Drop-in replacement for Snippy with feature parity.

    Examples:

        $ snippy-ng run --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.pipeline import Pipeline
    from snippy_ng.stages.setup import PrepareReference
    from snippy_ng.stages.clean_reads import CleanReadsFastp
    from snippy_ng.stages.alignment import BWAMEMReadsAligner, MinimapAligner, PreAlignedReads
    from snippy_ng.stages.alignment_filtering import AlignmentFilter
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
        
        # Clean reads (optional)
        if kwargs["clean_reads"] and kwargs["reads"]:
            clean_reads_stage = CleanReadsFastp(**kwargs)
            # Update reads to use cleaned reads
            kwargs["reads"] = [clean_reads_stage.output.cleaned_r1]
            if clean_reads_stage.output.cleaned_r2:
                kwargs["reads"].append(clean_reads_stage.output.cleaned_r2)
            stages.append(clean_reads_stage)
        
        # Aligner
        if kwargs["bam"]:
            aligner = PreAlignedReads(**kwargs)
        elif kwargs["aligner"] == "bwamem":
            aligner = BWAMEMReadsAligner(**kwargs)
        else:
            kwargs["aligner_opts"] = "-x sr " + kwargs.get("aligner_opts", "")
            aligner = MinimapAligner(**kwargs)
        kwargs["bam"] = aligner.output.bam
        stages.append(aligner)
        # Filter alignment
        align_filter = AlignmentFilter(**kwargs)
        kwargs["bam"] = align_filter.output.bam
        stages.append(align_filter)
        # SNP calling
        stages.append(FreebayesCaller(**kwargs))
    except ValidationError as e:
        error(e)
    
    # Move from CLI land into Pipeline land
    snippy = Pipeline(stages=stages)
    snippy.welcome()

    if not kwargs.get("skip_check", False):
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

