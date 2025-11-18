from pathlib import Path
import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@snippy_global_options
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--assembly", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Assembly in FASTA format")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
@click.option("--header", default=None, type=click.STRING, help="Header for the output FASTA file (if not provided, reference headers are kept)")
def asm(**config):
    """
    Assembly based SNP calling pipeline

    Examples:

        $ snippy-ng asm --reference ref.fa --assembly assembly.fa --outdir output
    """
    from snippy_ng.snippy import Snippy
    from snippy_ng.stages.filtering import VcfFilter
    from snippy_ng.exceptions import DependencyError, MissingOutputError
    from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
    from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
    from snippy_ng.stages.compression import BgzipCompressor
    from snippy_ng.stages.masks import ApplyMask, HetMask
    from snippy_ng.stages.copy import CopyFasta
    from snippy_ng.cli.utils import error
    from snippy_ng.cli.utils.common import load_or_prepare_reference
    from pydantic import ValidationError
    from snippy_ng.stages.alignment import AssemblyAligner
    from snippy_ng.stages.calling import PAFCaller


    # Choose stages to include in the pipeline
    try:
        stages = []
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=config["reference"],
            reference_prefix=config.get("prefix", "ref")
        )
        config["reference"] = setup.output.reference
        config["features"] = setup.output.gff
        config["reference_index"] = setup.output.reference_index
        stages.append(setup)
        
        # Aligner 
        aligner = AssemblyAligner(**config)
        stages.append(aligner)
        # Call variants
        caller = PAFCaller(
            paf=aligner.output.paf,
            ref_dict=setup.output.reference_dict,
            **config
        )
        stages.append(caller)
        # Filter VCF
        variant_filter = VcfFilter(
            vcf=caller.output.vcf,
            # hard code for asm-based calling
            min_depth=1,
            min_qual=60,
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
        # we should use config['fasta'] from now on
        config['reference'] = pseudo.output.fasta
        
        # Apply depth masking
        missing_mask = ApplyMask(
            fasta=config['reference'],
            mask_bed=caller.output.missing_bed,
            mask_char="-",
            **config
        )
        stages.append(missing_mask)
        config['reference'] = missing_mask.output.masked_fasta 

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
        snippy.run(quiet=config["quiet"], continue_last_run=config["continue_last_run"], keep_incomplete=config["keep_incomplete"])
    except MissingOutputError as e:
        snippy.error(e)
        return 1
    except RuntimeError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup()
    snippy.goodbye()
