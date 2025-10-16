from pathlib import Path
import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True})
@snippy_global_options
@click.option("--reference", "--ref", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--assembly", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Assembly in FASTA format")
@click.option("--aligner-opts", default='', type=click.STRING, help="Extra options for the aligner")
@click.option("--mask", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Mask file (BED format) to mask regions in the reference with Ns")
@click.option("--min-depth", default=1, type=click.INT, help="Minimum coverage to call a variant")
@click.option("--min-qual", default=60, type=click.FLOAT, help="Minimum QUAL threshold for heterozygous/low quality site masking")
@click.option("--prefix", default='snps', type=click.STRING, help="Prefix for output files")
@click.option("--header", default=None, type=click.STRING, help="Header for the output FASTA file (if not provided, reference headers are kept)")
def asm(**kwargs):
    """
    Assembly based SNP calling pipeline

    Examples:

        $ snippy-ng asm --reference ref.fa --assembly assembly.fa --outdir output
    """
    from snippy_ng.pipeline import Pipeline
    from snippy_ng.stages.filtering import VcfFilter
    from snippy_ng.exceptions import DependencyError, MissingOutputError
    from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
    from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
    from snippy_ng.stages.compression import BgzipCompressor
    from snippy_ng.stages.masks import ApplyMask, HetMask
    from snippy_ng.stages.copy import CopyFasta
    from snippy_ng.cli.utils import error
    from snippy_ng.cli.utils.reference import load_or_prepare_reference
    from pydantic import ValidationError
    from snippy_ng.stages.alignment import AssemblyAligner
    from snippy_ng.stages.calling import PAFCaller


    if not kwargs["force"] and kwargs["outdir"].exists():
        error(f"Output folder '{kwargs['outdir']}' already exists! Use --force to overwrite.")

    # check if output folder exists
    if not kwargs["outdir"].exists():
        kwargs["outdir"].mkdir(parents=True, exist_ok=True)

    # Choose stages to include in the pipeline
    try:
        stages = []
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=kwargs["reference"],
            reference_prefix=kwargs.get("prefix", "ref")
        )
        kwargs["reference"] = setup.output.reference
        kwargs["features"] = setup.output.gff
        kwargs["reference_index"] = setup.output.reference_index
        stages.append(setup)
        
        # Aligner 
        aligner = AssemblyAligner(**kwargs)
        stages.append(aligner)
        # Call variants
        caller = PAFCaller(
            paf=aligner.output.paf,
            ref_dict=setup.output.reference_dict,
            **kwargs
        )
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
        # we should use kwargs['fasta'] from now on
        kwargs['reference'] = pseudo.output.fasta
        
        # Apply depth masking
        missing_mask = ApplyMask(
            fasta=kwargs['reference'],
            mask_bed=caller.output.missing_bed,
            mask_char="-",
            **kwargs
        )
        stages.append(missing_mask)
        kwargs['reference'] = missing_mask.output.masked_fasta 

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