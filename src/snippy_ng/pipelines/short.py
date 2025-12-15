from pathlib import Path
from typing import Optional, List
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
from snippy_ng.pipelines.common import load_or_prepare_reference


def create_short_pipeline_stages(
    reference: str,
    reads: List[str],
    prefix: str = "snps",
    bam: Optional[str] = None,
    clean_reads: bool = False,
    downsample: Optional[float] = None,
    aligner: str = "minimap2",
    aligner_opts: str = "",
    freebayes_opts: str = "",
    mask: Optional[str] = None,
    min_depth: int = 10,
    min_qual: float = 100,
    tmpdir: Path = Path("/tmp"),
    cpus: int = 1,
    ram: int = 8,
) -> list:
    stages = []
    globals = {'prefix': prefix, 'cpus': cpus, 'ram': ram, 'tmpdir': tmpdir}
    
    # Setup reference (load existing or prepare new)
    setup = load_or_prepare_reference(
        reference_path=reference,
        reference_prefix=prefix,
    )
    reference_file = setup.output.reference
    features_file = setup.output.gff
    reference_index = setup.output.reference_index
    stages.append(setup)
    
    # Track current reads through potential cleaning and downsampling
    current_reads = reads.copy() if reads else []
    
    if downsample and current_reads:
        from snippy_ng.stages.downsample_reads import RasusaDownsampleReadsByCoverage
        from snippy_ng.at_run_time import get_genome_length
        
        # We need the genome length at run time (once we know the reference)
        genome_length = get_genome_length(setup.output.meta)
        downsample_stage = RasusaDownsampleReadsByCoverage(
            coverage=downsample,
            genome_length=genome_length,
            reads=current_reads,
            **globals
        )
        # Update reads to use downsampled reads
        current_reads = [downsample_stage.output.downsampled_r1]
        if downsample_stage.output.downsampled_r2:
            current_reads.append(downsample_stage.output.downsampled_r2)
        stages.append(downsample_stage)
    
    # Clean reads (optional)
    if clean_reads and current_reads:
        clean_reads_stage = FastpCleanReads(
            reads=current_reads,
            **globals
        )
        # Update reads to use cleaned reads
        current_reads = [clean_reads_stage.output.cleaned_r1]
        if clean_reads_stage.output.cleaned_r2:
            current_reads.append(clean_reads_stage.output.cleaned_r2)
        stages.append(clean_reads_stage)

    # Aligner
    if bam:
        aligner_stage = PreAlignedReads(
            bam=Path(bam),
            reference=reference_file,
            **globals
        )
    elif aligner == "bwamem":
        aligner_stage = BWAMEMReadsAligner(
            reads=current_reads,
            reference=reference_file,
            reference_index=reference_index,
            aligner_opts=aligner_opts,
            **globals
        )
    else:
        # Minimap2
        minimap_opts = "-x sr " + aligner_opts
        aligner_stage = MinimapAligner(
            reads=current_reads,
            reference=reference_file,
            aligner_opts=minimap_opts,
            **globals
        )
    
    if not bam:
        # SeqKit read statistics
        stats_stage = SeqKitReadStatsBasic(
            reads=current_reads,
            **globals
        )
        stages.append(stats_stage)
    
    current_bam = aligner_stage.output.bam
    stages.append(aligner_stage)
    
    # Filter alignment
    align_filter = BamFilter(
        bam=current_bam,
        **globals
    )
    current_bam = align_filter.output.bam
    stages.append(align_filter)
    
    # SNP calling
    caller = FreebayesCaller(
        bam=current_bam,
        reference=reference_file,
        reference_index=reference_index,
        fbopt=freebayes_opts,
        **globals
    )
    stages.append(caller)
    
    # Filter VCF
    variant_filter = VcfFilter(
        vcf=caller.output.vcf,
        reference=reference_file,
        min_depth=min_depth,
        min_qual=min_qual,
        **globals
    )
    stages.append(variant_filter)
    variants_file = variant_filter.output.vcf
    
    # Consequences calling
    consequences = BcftoolsConsequencesCaller(
        variants=variants_file,
        features=features_file,
        reference=reference_file,
        **globals
    )
    stages.append(consequences)
    
    # Compress VCF
    gzip = BgzipCompressor(
        input=consequences.output.annotated_vcf,
        suffix="gz",
        **globals
    )
    stages.append(gzip)
    
    # Pseudo-alignment
    pseudo = BcftoolsPseudoAlignment(
        vcf_gz=gzip.output.compressed,
        reference=reference_file,
        **globals
    )
    stages.append(pseudo)
    
    # Track the current reference/fasta through the masking stages
    current_fasta = pseudo.output.fasta
    
    # Apply depth masking
    depth_mask = DepthMask(
        bam=current_bam,
        fasta=current_fasta,
        min_depth=min_depth,
        **globals
    )
    stages.append(depth_mask)
    current_fasta = depth_mask.output.masked_fasta

    # Apply heterozygous and low quality sites masking
    het_mask = HetMask(
        vcf=caller.output.vcf,  # Use raw VCF for complete site information
        fasta=current_fasta,
        min_qual=min_qual,
        **globals
    )
    stages.append(het_mask)
    current_fasta = het_mask.output.masked_fasta
    
    # Apply user mask if provided
    if mask:
        user_mask = ApplyMask(
            fasta=current_fasta,
            mask_bed=Path(mask),
            **globals
        )
        stages.append(user_mask)
        current_fasta = user_mask.output.masked_fasta

    # Copy final masked consensus to standard output location
    copy_final = CopyFasta(
        input=current_fasta,
        output_path=f"{prefix}.pseudo.fna",
        **globals
    )
    stages.append(copy_final)

    return stages
    


