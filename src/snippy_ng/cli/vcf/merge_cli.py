from pathlib import Path

import click

from snippy_ng.utils.vcf_merge import VcfMergeError, merge_vcfs, read_vcf_file_list, write_upset_svg


@click.command("merge", context_settings={"show_default": True})
@click.argument("vcfs", nargs=-1, type=click.Path(path_type=Path, exists=True, readable=True, dir_okay=False))
@click.option("--file-list", "file_lists", multiple=True, type=click.Path(path_type=Path, exists=True, readable=True, dir_okay=False), help="File with one VCF path per line; may be supplied more than once")
@click.option("--output", "-o", "output_vcf", type=click.Path(path_type=Path, dir_okay=False), help="Merged VCF/VCF.GZ/BCF output (defaults to stdout)")
@click.option("--upset-plot", type=click.Path(path_type=Path, dir_okay=False), help="Optional SVG UpSet plot of input VCF overlap")
@click.option("--max-intersections", type=click.IntRange(min=1), default=40, help="Maximum intersections to show in the UpSet plot")
@click.option("--force-samples/--no-force-samples", default=True, help="Keep duplicate sample names as separate columns")
def merge(
    vcfs: tuple[Path, ...],
    file_lists: tuple[Path, ...],
    output_vcf: Path | None,
    upset_plot: Path | None,
    max_intersections: int,
    force_samples: bool,
) -> None:
    """Merge coordinate-sorted VCFs with ``bcftools merge --no-index``."""
    input_vcfs = [*vcfs, *read_vcf_file_list(file_lists)]
    missing = [path for path in input_vcfs if not path.is_file()]
    if missing:
        raise click.UsageError(f"VCF does not exist: {missing[0]}")
    if len(input_vcfs) < 2:
        raise click.UsageError("Provide at least two VCFs, directly or with --file-list")
    try:
        merge_vcfs(input_vcfs, output_vcf, force_samples=force_samples)
    except VcfMergeError as exc:
        raise click.ClickException(str(exc)) from exc
    if upset_plot is not None:
        write_upset_svg(input_vcfs, upset_plot, max_intersections=max_intersections)
