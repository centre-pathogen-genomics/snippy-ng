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
@click.option("--freebayes-opts", default='', type=click.STRING, help="Extra options for Freebayes")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--min-depth", default=10, type=click.INT, help="Minimum coverage to call a variant")
@click.option("--min-qual", default=100, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
@click.option("--header", default=None, type=click.STRING, help="Header for the output FASTA file (if not provided, reference headers are kept)")
def short(
    *,
    # Global options (from snippy_global_options decorator)
    outdir: Path,
    tmpdir: Path,
    cpus: int,
    ram: int,
    force: bool,
    quiet: bool,
    debug: bool,
    skip_check: bool,
    check: bool,
    keep_incomplete: bool,
    continue_last_run: bool,
    # Short command specific options
    reference: Path,
    r1: Path | None,
    r2: Path | None,
    bam: Path | None,
    clean_reads: bool,
    downsample: float | None,
    aligner: str,
    aligner_opts: str,
    freebayes_opts: str,
    mask: Path | None,
    min_depth: int,
    min_qual: float,
    prefix: str,
    header: str | None,
):
    """
    Short read based SNP calling pipeline

    Examples:

        $ snippy-ng short --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.snippy import Snippy
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
    from snippy_ng.cli.utils.common import load_or_prepare_reference
    from pydantic import ValidationError
    
    # Collect all arguments into a dictionary
    config = {
        "outdir": outdir,
        "tmpdir": tmpdir,
        "cpus": cpus,
        "ram": ram,
        "force": force,
        "quiet": quiet,
        "debug": debug,
        "skip_check": skip_check,
        "check": check,
        "continue": continue_last_run,
        "keep_incomplete": keep_incomplete,
        "reference": reference,
        "r1": r1,
        "r2": r2,
        "bam": bam,
        "clean_reads": clean_reads,
        "downsample": downsample,
        "aligner": aligner,
        "aligner_opts": aligner_opts,
        "freebayes_opts": freebayes_opts,
        "mask": mask,
        "min_depth": min_depth,
        "min_qual": min_qual,
        "prefix": prefix,
        "header": header,
    }

    # combine R1 and R2 into reads
    config["reads"] = []
    


    # combine R1 and R2 into reads
    config["reads"] = []
    if config.get("r1"):
        config["reads"].append(config["r1"])
    if config.get("r2"):
        config["reads"].append(config["r2"])
    if not config["reads"] and not config.get("bam"):
        error("Please provide reads or a BAM file!")
    
    
    # Choose stages to include in the pipeline
    try:
        stages = []
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=config["reference"],
            reference_prefix=config.get("prefix", "ref"),
        )
        config["reference"] = setup.output.reference
        config["features"] = setup.output.gff
        config["reference_index"] = setup.output.reference_index
        stages.append(setup)
        
        # Clean reads (optional)
        if config["clean_reads"] and config["reads"]:
            clean_reads_stage = FastpCleanReads(**config)
            # Update reads to use cleaned reads
            config["reads"] = [clean_reads_stage.output.cleaned_r1]
            if clean_reads_stage.output.cleaned_r2:
                config["reads"].append(clean_reads_stage.output.cleaned_r2)
            stages.append(clean_reads_stage)
        if config.get("downsample") and config.get("reads"):
            from snippy_ng.stages.downsample_reads import RasusaDownsampleReadsByCoverage
            from snippy_ng.runtime import at_run_time, genome_length_getter
            
            # We need the genome length at run time (once we know the reference)
            genome_length=at_run_time(genome_length_getter(setup.output.meta))
            downsample_stage = RasusaDownsampleReadsByCoverage(
                coverage=config["downsample"],
                genome_length=genome_length,
                **config
            )
            # Update reads to use downsampled reads
            config["reads"] = [downsample_stage.output.downsampled_r1]
            if downsample_stage.output.downsampled_r2:
                config["reads"].append(downsample_stage.output.downsampled_r2)
            stages.append(downsample_stage)
        
        # Aligner
        if config["bam"]:
            aligner = PreAlignedReads(**config)
        elif config["aligner"] == "bwamem":
            aligner = BWAMEMReadsAligner(**config)
        else:
            config["aligner_opts"] = "-x sr " + config.get("aligner_opts", "")
            aligner = MinimapAligner(**config)
        if not config["bam"]:
            # SeqKit read statistics
            stages.append(SeqKitReadStatsBasic(**config))
        config["bam"] = aligner.output.bam
        stages.append(aligner)
        # Filter alignment
        align_filter = BamFilter(**config)
        config["bam"] = align_filter.output.bam
        stages.append(align_filter)
        # SNP calling
        caller = FreebayesCaller(
            fbopt=config["freebayes_opts"],
            **config
        )
        stages.append(caller)
        # Filter VCF
        variant_filter = VcfFilter(
            vcf=caller.output.vcf,
            **config,
        )
        stages.append(variant_filter)
        config["variants"] = variant_filter.output.vcf
        # Consequences calling
        consequences = BcftoolsConsequencesCaller(**config) 
        stages.append(consequences)
        # Compress VCF
        gzip = BgzipCompressor(
            input=consequences.output.annotated_vcf,
            suffix="gz",
            **config,
        )
        stages.append(gzip)
        # Pseudo-alignment
        pseudo = BcftoolsPseudoAlignment(vcf_gz=gzip.output.compressed, **config)
        stages.append(pseudo)
        config['reference'] = pseudo.output.fasta
        
        # Apply depth masking
        depth_mask = DepthMask(
            **config
        )
        stages.append(depth_mask)
        config['reference'] = depth_mask.output.masked_fasta

        # Apply heterozygous and low quality sites masking
        het_mask = HetMask(
            vcf=caller.output.vcf,  # Use raw VCF for complete site information
            **config
        )
        stages.append(het_mask)
        config['reference'] = het_mask.output.masked_fasta
        
        # Apply user mask if provided
        if config["mask"]:
            user_mask = ApplyMask(
                mask_bed=Path(config["mask"]),
                **config
            )
            stages.append(user_mask)
            config['reference'] = user_mask.output.masked_fasta

        # Copy final masked consensus to standard output location
        from snippy_ng.stages.copy import CopyFasta
        copy_final = CopyFasta(
            input=config['reference'],
            output_path=f"{config['prefix']}.pseudo.fna",
            **config,
        )
        stages.append(copy_final)
            
    except ValidationError as e:
        error(e)
    
    # Move from CLI land into Pipeline land
    snippy = Snippy(stages=stages)
    snippy.welcome()

    if not config.get("skip_check", False):
        try:
            snippy.validate_dependencies()
        except DependencyError as e:
            snippy.error(f"Invalid dependencies! Please install '{e}' or use --skip-check to ignore.")
            return 1
    
    if config["check"]:
        return 0

    # Set working directory to output folder
    snippy.set_working_directory(config["outdir"])
    try:
        snippy.run(quiet=config["quiet"], continue_last_run=config["continue"], keep_incomplete=config["keep_incomplete"])
    except MissingOutputError as e:
        snippy.error(e)
        return 1
    except RuntimeError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup()
    snippy.goodbye()


