from pathlib import Path

import click

from snippy_ng.cli.utils import AbsolutePath


@click.command(context_settings={"show_default": True})
@click.argument("alignment", required=True, type=AbsolutePath(exists=True, readable=True))
@click.option("--gff", type=AbsolutePath(exists=True, readable=True), help="Reference GFF for per-feature CNV estimates")
@click.option("--feature", "feature_type", default="gene", help="GFF feature type to use when --gff is supplied")
@click.option("--known-single-copy", help="Single-copy baseline region as START,END on the largest contig or CONTIG:START-END")
@click.option("--no-header", is_flag=True, default=False, help="Do not print the TSV header")
def cnv(
    alignment: Path,
    gff: Path | None,
    feature_type: str,
    known_single_copy: str | None,
    no_header: bool,
):
    """
    Estimate copy number variation from aligned BAM or CRAM depth.

    The largest contig is assumed to have one copy number. Other contigs are assigned
    rounded integer copy numbers from their depth relative to that baseline.
    """
    import subprocess
    import sys

    from snippy_ng.utils.cnv import (
        CNVError,
        copy_number_variation,
        write_cnv_table,
        write_feature_cnv_table,
    )

    try:
        rows = copy_number_variation(
            alignment=alignment,
            gff=gff,
            feature_type=feature_type,
            known_single_copy=known_single_copy,
        )
        if gff:
            write_feature_cnv_table(rows, sys.stdout, include_header=not no_header)
        else:
            write_cnv_table(rows, sys.stdout, include_header=not no_header)
    except subprocess.CalledProcessError as exc:
        message = exc.stderr.strip() or str(exc)
        raise click.ClickException(f"samtools failed: {message}") from exc
    except CNVError as exc:
        raise click.ClickException(str(exc)) from exc
