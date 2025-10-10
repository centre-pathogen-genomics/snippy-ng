from pathlib import Path
import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True}, short_help="Run SNP calling pipeline for short reads")
@snippy_global_options
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--R1", "--pe1", "--left", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R2 (right)")
@click.option("--clean-reads", is_flag=True, default=False, help="Clean and filter reads with fastp before alignment")
@click.option("--downsample", type=click.FLOAT, default=None, help="Downsample reads to a specified coverage (e.g., 30.0 for 30x coverage)")
@click.option("--aligner", default="minimap2", type=click.Choice(["minimap2", "bwamem"]), help="Aligner program to use")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--min-depth", default=10, type=click.INT, help="Minimum coverage to call a variant")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
def short(**kwargs):
    """
    Drop-in replacement for Snippy with feature parity.

    Examples:

        $ snippy-ng short --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output
    """
    from snippy_ng.pipeline import Pipeline
    from snippy_ng.stages.setup import LoadReference, PrepareReference
    from snippy_ng.stages.clean_reads import FastpCleanReads
    from snippy_ng.stages.stats import SeqKitReadStatsBasic
    from snippy_ng.stages.alignment import BWAMEMReadsAligner, MinimapAligner, PreAlignedReads
    from snippy_ng.stages.alignment_filtering import AlignmentFilter
    from snippy_ng.stages.calling import FreebayesCaller
    from snippy_ng.exceptions import DependencyError, MissingOutputError
    from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
    from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
    from snippy_ng.stages.compression import BgzipCompressor
    from snippy_ng.stages.copy import CopyFile
    from snippy_ng.seq_utils import guess_format
    from snippy_ng.cli.utils import error
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
    stages = []
    try:
        if Path(kwargs["reference"]).is_dir():
            setup = LoadReference(
                    reference_dir=kwargs["reference"], 
                    **kwargs,
                )
            kwargs["reference"] = setup.output.reference
            kwargs["features"] = setup.output.gff
            kwargs["reference_index"] = setup.output.reference_index
            stages.append(setup)
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
            
            # Create a closure that captures the setup stage
            def make_genome_length_getter():
                _setup = setup if 'setup' in locals() else None
                _outdir = kwargs["outdir"]
                
                def get_genome_length():
                    import json
                    # Use the setup stage's metadata file if available
                    if _setup and hasattr(_setup, 'output'):
                        meta_path = _setup.output.meta
                    else:
                        meta_path = _outdir / "reference" / "ref.json"
                    with open(meta_path, 'r') as f:
                        metadata = json.load(f)
                    return int(metadata['total_length'])
                
                return get_genome_length
            
            downsample_stage = RasusaDownsampleReadsByCoverage(
                coverage=kwargs["downsample"],
                genome_length=at_run_time(make_genome_length_getter()),
                **kwargs
            )
            # Update reads to use downsampled reads
            kwargs["reads"] = [downsample_stage.output.downsampled_r1]
            if downsample_stage.output.downsampled_r2:
                kwargs["reads"].append(downsample_stage.output.downsampled_r2)
            stages.append(downsample_stage)
        # SeqKit read statistics
        stages.append(SeqKitReadStatsBasic(**kwargs))
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
        caller = FreebayesCaller(**kwargs)
        stages.append(caller)
        kwargs["variants"] = caller.output.filter_vcf
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
        kwargs["reference"] = pseudo.output.fasta
        
        # Apply masking stages
        masking_stages = apply_masks(kwargs)
        stages.extend(masking_stages)

        # rename final consensus output to standard prefix
        final_consensus = masking_stages[-1].output.fasta if masking_stages else pseudo.output.fasta
        copy_final = CopyFile(
            input=final_consensus,
            output_path=f"{kwargs['prefix']}.fasta",
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
        snippy.run(quiet=kwargs["quiet"])
    except MissingOutputError as e:
        snippy.error(e)
        return 1
    except RuntimeError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup()
    snippy.goodbye()


def apply_masks(kwargs):
    """
    Apply masking stages in sequence: min-depth mask (optional), zero-depth mask, and user mask (optional).
    
    Args:
        kwargs: Dictionary containing pipeline parameters
        stages: List of existing pipeline stages (for reference)
    
    Returns:
        List of masking stages to add to the pipeline
    """
    from snippy_ng.stages.masks import ZeroDepthMask, ApplyMask, MinDepthMask
    
    masking_stages = []
    original_prefix = kwargs["prefix"]
    
    # Define masking configurations
    mask_configs = []
    
    # Min depth mask (conditional)
    if kwargs["min_depth"] > 0:
        mask_configs.append({
            "mask_stage": MinDepthMask,
            "mask_kwargs": {"filter": f"<{kwargs['min_depth']}"},
            "apply_kwargs": {"char": "n"},
            "suffix": "mindepth"
        })
    
    # Zero depth mask (always applied)
    mask_configs.append({
        "mask_stage": ZeroDepthMask,
        "mask_kwargs": {},
        "apply_kwargs": {"char": "-"},
        "suffix": "zero_depth"
    })
    
    # User mask (conditional)
    if kwargs["mask"]:
        mask_configs.append({
            "mask_stage": None,  # No mask generation stage needed
            "bed_file": Path(kwargs["mask"]),
            "apply_kwargs": {"char": "N"},
            "suffix": "user"
        })
    
    # Apply masks in sequence
    for config in mask_configs:
        # Generate mask if needed
        if config["mask_stage"]:
            mask_stage_kwargs = {**kwargs, **config.get("mask_kwargs", {})}
            mask_stage = config["mask_stage"](**mask_stage_kwargs)
            masking_stages.append(mask_stage)
            bed_file = mask_stage.output.bed
        else:
            bed_file = config["bed_file"]
        
        # Apply mask to reference
        kwargs["prefix"] = f"{original_prefix}.{config['suffix']}"
        apply_mask_kwargs = {**kwargs, "bed": bed_file, **config["apply_kwargs"]}
        apply_mask = ApplyMask(**apply_mask_kwargs)
        masking_stages.append(apply_mask)
        
        # Update reference for next iteration
        kwargs["reference"] = apply_mask.output.fasta
        kwargs["prefix"] = original_prefix
    
    return masking_stages
