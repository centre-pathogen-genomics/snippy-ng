from pathlib import Path

import click

from snippy_ng.cli.utils import AbsolutePath



@click.command("strand-bias", context_settings={"show_default": True})
@click.option("--alignment", "--bam", "--cram", required=True, type=AbsolutePath(exists=True, readable=True), help="Input BAM or CRAM used to derive strand counts")
@click.option("--reference", "--ref", required=True, type=AbsolutePath(exists=True, readable=True), help="Reference FASTA for samtools mpileup")
@click.option("--output", "-o", required=False, type=click.Path(path_type=Path, dir_okay=False, writable=True), help="Output VCF path. Defaults to stdout")
@click.option("--filter-pvalue", type=click.FloatRange(min=0.0, max=1.0), default=None, help="Optionally mark variants with FILTER=StrandBias when the Fisher p-value is below this threshold")
@click.argument("vcf", type=AbsolutePath(exists=True, readable=True))
def strand_bias(alignment: Path, reference: Path, output: Path | None, filter_pvalue: float | None, vcf: Path):
    """Annotate VCF records with strand-bias counts and Fisher exact p-values."""
    import subprocess

    from snippy_ng.utils.strand_bias import StrandBiasError, annotate_vcf_strand_bias

    try:
        annotate_vcf_strand_bias(
            vcf=vcf,
            bam=alignment,
            reference=reference,
            output=output.absolute() if output is not None else None,
            filter_pvalue=filter_pvalue,
        )
    except subprocess.CalledProcessError as exc:
        message = exc.stderr.strip() or str(exc)
        raise click.ClickException(f"samtools failed: {message}") from exc
    except StrandBiasError as exc:
        raise click.ClickException(str(exc)) from exc
