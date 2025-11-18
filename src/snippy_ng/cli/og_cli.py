import click
from snippy_ng.cli.utils.globals import CommandWithGlobals, snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={'show_default': True}, short_help="Backwards-compatible SNP calling pipeline", deprecated=True)
@click.option("--reference", required=True, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reference genome (FASTA or GenBank)")
@click.option("--R1", "--pe1", "--left", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default=None, type=click.Path(exists=True, resolve_path=True, readable=True), help="Reads, paired-end R2 (right)")
@click.option("--se", "--single", default='', type=click.STRING, help="Single-end reads")
@click.option("--ctgs", "--contigs", default='', type=click.STRING, help="Use these contigs instead of reads")
@click.option("--peil", default='', type=click.STRING, help="Paired-end interleaved reads")
@click.option("--bam", default=None, type=click.Path(exists=True, resolve_path=True), help="Use this BAM file instead of aligning reads")
@click.option("--targets", default='', type=click.STRING, help="Only call SNPs from this BED file")
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
    Drop-in replacement for Snippy that maps command-line arguments to the new Snippy-NG short pipeline.

    Note that although the interface is the same as Snippy it will produce slightly different results
    due to differences in the underlying tools and versions used.
    
    Examples:\n
        $ alias snippy='snippy-ng og'\n
        $ snippy --reference ref.fa --R1 reads_1.fq --R2 reads_2.fq --outdir output

    SYNOPSIS
      snippy 4.6.0 - fast bacterial variant calling from NGS reads
    USAGE
      snippy [options] --outdir <dir> --ref <ref> --R1 <R1.fq.gz> --R2 <R2.fq.gz>
      snippy [options] --outdir <dir> --ref <ref> --ctgs <contigs.fa>
      snippy [options] --outdir <dir> --ref <ref> --bam <reads.bam>
    GENERAL
      --help           This help
      --version        Print version and exit
      --citation       Print citation for referencing snippy
      --check          Check dependences are installed then exit (default OFF)
      --force          Force overwrite of existing output folder (default OFF)
      --quiet          No screen output (default OFF)
    RESOURCES
      --cpus N         Maximum number of CPU cores to use (default '8')
      --ram N          Try and keep RAM under this many GB (default '8')
      --tmpdir F       Fast temporary storage eg. local SSD (default '/var/folders/hs/3sl81nqd6mzcbz1sc_td3bv00000gn/T/')
    INPUT
      --reference F    Reference genome. Supports FASTA, GenBank, EMBL (not GFF) (default '')
      --R1 F           Reads, paired-end R1 (left) (default '')
      --R2 F           Reads, paired-end R2 (right) (default '')
      --se F           Single-end reads (default '')
      --ctgs F         Don't have reads use these contigs (default '')
      --peil F         Reads, paired-end R1/R2 interleaved (default '')
      --bam F          Use this BAM file instead of aligning reads (default '')
      --targets F      Only call SNPs from this BED file (default '')
      --subsample n.n  Subsample FASTQ to this proportion (default '1')
    OUTPUT
      --outdir F       Output folder (default '')
      --prefix F       Prefix for output files (default 'snps')
      --report         Produce report with visual alignment per variant (default OFF)
      --cleanup        Remove most files not needed for snippy-core (inc. BAMs!) (default OFF)
      --rgid F         Use this @RG ID: in the BAM header (default '')
      --unmapped       Keep unmapped reads in BAM and write FASTQ (default OFF)
    PARAMETERS
      --mapqual N      Minimum read mapping quality to consider (default '60')
      --basequal N     Minimum base quality to consider (default '13')
      --mincov N       Minimum site depth to for calling alleles (default '10')
      --minfrac n.n    Minimum proportion for variant evidence (0=AUTO) (default '0')
      --minqual n.n    Minimum QUALITY in VCF column 6 (default '100')
      --maxsoft N      Maximum soft clipping to allow (default '10')
      --bwaopt F       Extra BWA MEM options, eg. -x pacbio (default '')
      --fbopt F        Extra Freebayes options, eg. --theta 1E-6 --read-snp-limit 2 (default '')
    """
    from snippy_ng.cli.short_cli import short

    short.callback(
        reference=kwargs["reference"],
        r1=kwargs["r1"],
        r2=kwargs["r2"],
        bam=kwargs["bam"],
        prefix=kwargs["prefix"],
        min_depth=kwargs["mincov"],
        min_qual=kwargs["minqual"],
        mask=kwargs["targets"],
        header=None,
        aligner="bwamem",
        aligner_opts=kwargs["bwaopt"], 
    )

