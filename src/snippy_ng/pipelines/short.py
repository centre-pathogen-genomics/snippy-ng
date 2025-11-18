from pathlib import Path
from snippy_ng.stages.clean_reads import FastpCleanReads
from snippy_ng.stages.stats import SeqKitReadStatsBasic
from snippy_ng.stages.alignment import BWAMEMReadsAligner, MinimapAligner, PreAlignedReads
from snippy_ng.stages.filtering import BamFilter, VcfFilter
from snippy_ng.stages.calling import FreebayesCaller
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import BgzipCompressor
from snippy_ng.stages.masks import DepthMask, ApplyMask, HetMask
from snippy_ng.stages.copy import CopyFasta
from snippy_ng.cli.utils.common import load_or_prepare_reference


def create_short_stages(config: dict) -> list:
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
        from snippy_ng.at_run_time import get_genome_length
        
        # We need the genome length at run time (once we know the reference)
        genome_length=get_genome_length(setup.output.meta)
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
    copy_final = CopyFasta(
        input=config['reference'],
        output_path=f"{config['prefix']}.pseudo.fna",
        **config,
    )
    stages.append(copy_final)

    return stages
    


