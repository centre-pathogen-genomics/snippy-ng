from pathlib import Path
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
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--min-depth", default=10, type=click.INT, help="Minimum coverage to call a variant")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
@click.option("--header", default=None, type=click.STRING, help="Header for the output FASTA file (if not provided, reference headers are kept)")
def short(**kwargs):
    """
    Short read based SNP calling pipeline

    Examples:

        $ snippy-ng short --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.pipeline import Pipeline
    from snippy_ng.stages.clean_reads import FastpCleanReads
    from snippy_ng.stages.stats import SeqKitReadStatsBasic
    from snippy_ng.stages.alignment import BWAMEMReadsAligner, MinimapAligner, PreAlignedReads
    from snippy_ng.stages.filtering import BamFilter, VcfFilter
    from snippy_ng.stages.calling import FreebayesCaller
    from snippy_ng.exceptions import DependencyError, MissingOutputError
    from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
    from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
    from snippy_ng.stages.compression import BgzipCompressor
    from snippy_ng.stages.masks import DepthMask, ApplyMask, HetMask
    from snippy_ng.stages.copy import CopyFasta
    from snippy_ng.cli.utils import error
    from snippy_ng.cli.utils.reference import load_or_prepare_reference
    from pydantic import ValidationError


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
    try:
        stages = []
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=kwargs["reference"],
            reference_prefix=kwargs.get("prefix", "ref"),
        )
        kwargs["reference"] = setup.output.reference
        kwargs["features"] = setup.output.gff
        kwargs["reference_index"] = setup.output.reference_index
        stages.append(setup)
        
        # Clean reads (optional)
        if kwargs["clean_reads"] and kwargs["reads"]:
            clean_reads_stage = FastpCleanReads(**kwargs)
            # Update reads to use cleaned reads
            kwargs["reads"] = [clean_reads_stage.output.cleaned_r1]
            if clean_reads_stage.output.cleaned_r2:
                kwargs["reads"].append(clean_reads_stage.output.cleaned_r2)
            stages.append(clean_reads_stage)
        if kwargs.get("downsample"):
            from snippy_ng.stages.downsample_reads import RasusaDownsampleReadsByCoverage
            from snippy_ng.stages import at_run_time
            
            # We need the genome length at run time (once we know the reference)
            genome_length=at_run_time(genome_length_getter(setup.output.meta))
            downsample_stage = RasusaDownsampleReadsByCoverage(
                coverage=kwargs["downsample"],
                genome_length=genome_length,
                **kwargs
            )
            # Update reads to use downsampled reads
            kwargs["reads"] = [downsample_stage.output.downsampled_r1]
            if downsample_stage.output.downsampled_r2:
                kwargs["reads"].append(downsample_stage.output.downsampled_r2)
            stages.append(downsample_stage)
        
        # Aligner
        if kwargs["bam"]:
            aligner = PreAlignedReads(**kwargs)
        elif kwargs["aligner"] == "bwamem":
            aligner = BWAMEMReadsAligner(**kwargs)
        else:
            kwargs["aligner_opts"] = "-x sr " + kwargs.get("aligner_opts", "")
            aligner = MinimapAligner(**kwargs)
        if not kwargs["bam"]:
            # SeqKit read statistics
            stages.append(SeqKitReadStatsBasic(**kwargs))
        kwargs["bam"] = aligner.output.bam
        stages.append(aligner)
        # Filter alignment
        align_filter = BamFilter(**kwargs)
        kwargs["bam"] = align_filter.output.bam
        stages.append(align_filter)
        # SNP calling
        caller = FreebayesCaller(**kwargs)
        stages.append(caller)
        # Filter VCF
        variant_filter = VcfFilter(
            vcf=caller.output.vcf,
            **kwargs,
        )
        stages.append(variant_filter)
        kwargs["variants"] = variant_filter.output.vcf
        # Consequences calling
        consequences = BcftoolsConsequencesCaller(**kwargs) 
        stages.append(consequences)
        # Compress VCF
        gzip = BgzipCompressor(
            input=consequences.output.annotated_vcf,
            suffix="gz",
            **kwargs,
        )
        stages.append(gzip)
        # Pseudo-alignment
        pseudo = BcftoolsPseudoAlignment(vcf_gz=gzip.output.compressed, **kwargs)
        stages.append(pseudo)
        kwargs['reference'] = pseudo.output.fasta
        
        # Apply depth masking
        depth_mask = DepthMask(
            **kwargs
        )
        stages.append(depth_mask)
        kwargs['reference'] = depth_mask.output.masked_fasta

        # Apply heterozygous and low quality sites masking
        het_mask = HetMask(
            vcf=caller.output.vcf,  # Use raw VCF for complete site information
            **kwargs
        )
        stages.append(het_mask)
        kwargs['reference'] = het_mask.output.masked_fasta
        
        # Apply user mask if provided
        if kwargs["mask"]:
            user_mask = ApplyMask(
                mask_bed=Path(kwargs["mask"]),
                **kwargs
            )
            stages.append(user_mask)
            kwargs['reference'] = user_mask.output.masked_fasta

        # Copy final masked consensus to standard output location
        from snippy_ng.stages.copy import CopyFasta
        copy_final = CopyFasta(
            input=kwargs['reference'],
            output_path=f"{kwargs['prefix']}.pseudo.fna",
            **kwargs,
        )
        stages.append(copy_final)
            
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
        snippy.run(quiet=kwargs["quiet"], continue_last_run=kwargs["continue"], keep_incomplete=kwargs["keep_incomplete"])
    except MissingOutputError as e:
        snippy.error(e)
        return 1
    except RuntimeError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup()
    snippy.goodbye()


def genome_length_getter(reference_metadata: Path):
    """
    Because we don't know the genome length until run time (it depends on the reference provided),
    we create a closure that captures the setup stage and output directory, and returns a function
    that reads the genome length from the metadata file at run time.
    """
    def wraps():
        import json
        # Use the setup stage's metadata file if available
        with open(reference_metadata, 'r') as f:
            metadata = json.load(f)
        return int(metadata['total_length'])
    
    return wraps