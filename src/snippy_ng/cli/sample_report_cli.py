from pathlib import Path
from typing import Any, Optional

import click

from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import (
    CommandWithGlobals,
    GlobalOption,
    add_snippy_global_options,
    check_outdir_callback,
)


@click.command("sample-report", cls=CommandWithGlobals, context_settings={"show_default": True})
@click.option("--outdir", "-o", default="report", required=False, type=click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=Path), help="Output directory for the sample report", callback=check_outdir_callback, cls=GlobalOption)
@click.option("--prefix", "-p", default="sample-report", help="Prefix for HTML report", cls=GlobalOption)
@add_snippy_global_options(exclude=["prefix", "outdir"])
@click.argument("vcf", required=True, type=AbsolutePath(exists=True, readable=True))
@click.option("--alignment", "--cram", "--bam", required=False, type=AbsolutePath(exists=True, readable=True), help="Optional BAM or CRAM alignment file to embed after windowing")
@click.option("--reference", "--ref", required=False, type=AbsolutePath(exists=True, readable=True), help="Reference FASTA used by the alignment")
@click.option("--variant-scope", default="pass", type=click.Choice(["pass", "all"]), help="Variants to include in the report")
@click.option("--window-size", default=100, type=click.IntRange(min=0), help="Base pairs of context around each variant")
@click.option("--title", required=False, type=click.STRING, default="Snippy-NG Sample Report", help="Title for the HTML report")
@click.option("--sample-name", required=False, type=click.STRING, help="Optional sample name override")
def sample_report(
    vcf: Path,
    alignment: Path,
    reference: Path,
    variant_scope: str,
    window_size: int,
    title: str,
    sample_name: Optional[str],
    outdir: Optional[Path],
    prefix: str,
    **context: Any,
):
    """Create an HTML report with a variant table and optional embedded IGV.js alignment view."""
    from snippy_ng.context import Context
    from snippy_ng.pipelines.reports import SampleReportPipelineBuilder

    if alignment and not reference:
        raise click.UsageError("--reference is required when --alignment is provided")

    pipeline = SampleReportPipelineBuilder(
        vcf=vcf,
        alignment=alignment,
        reference=reference,
        title=title,
        sample_name=sample_name,
        variant_scope=variant_scope,
        window_size=window_size,
        prefix=prefix,
    ).build()

    context["outdir"] = outdir or Path(".").absolute()
    run_ctx = Context(**context)
    return pipeline.run(run_ctx)
