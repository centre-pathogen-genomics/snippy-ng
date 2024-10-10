import click

from snippy_ng.__about__ import __version__, EXE, AUTHOR, URL

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

def usage(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"Usage: {EXE} [OPTIONS]")
    ctx.exit()

@click.group(context_settings={"help_option_names": ["-h", "--help"]}, invoke_without_command=True)
@click.version_option(version=__version__, prog_name=EXE)
def snippy_ng():
    pass

@snippy_ng.command()
@click.option("--citation", is_flag=True, callback=show_citation, expose_value=False, help="Print citation for referencing this tool")
@click.option("--check/--no-check", default=False, help="Check dependencies are installed then exit")
@click.option("--force/--no-force", default=False, help="Force overwrite of existing output folder")
@click.option("--quiet/--no-quiet", default=False, help="No screen output")
@click.option("--cpus", default=8, type=int, help="Maximum number of CPU cores to use")
@click.option("--ram", default=8, type=int, help="Try and keep RAM under this many GB")
@click.option("--tmpdir", default='/tmp', type=click.Path(), help="Fast temporary storage eg. local SSD")
@click.option("--reference", default='', type=click.STRING, help="Reference genome (FASTA, GenBank, EMBL)")
@click.option("--R1", "--pe1", "--left", default='', type=click.STRING, help="Reads, paired-end R1 (left)")
@click.option("--R2", "--pe2", "--right", default='', type=click.STRING, help="Reads, paired-end R2 (right)")
@click.option("--se", "--single", default='', type=click.STRING, help="Single-end reads")
@click.option("--ctgs", "--contigs", default='', type=click.STRING, help="Use these contigs instead of reads")
@click.option("--peil", default='', type=click.STRING, help="Paired-end interleaved reads")
@click.option("--bam", default='', type=click.STRING, help="Use this BAM file instead of aligning reads")
@click.option("--targets", default='', type=click.STRING, help="Only call SNPs from this BED file")
@click.option("--subsample", default=1.0, type=float, help="Subsample FASTQ to this proportion")
@click.option("--outdir", default='', type=click.STRING, help="Output folder")
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
def snippy(**kwargs):
    """Snippy-NG: A tool for calling SNPs and other genomic variants."""
    click.echo("Running snippy-ng with options:")
    click.echo(kwargs)
    # Implement the core functionality based on the options

if __name__ == "__main__":
    snippy_ng()
