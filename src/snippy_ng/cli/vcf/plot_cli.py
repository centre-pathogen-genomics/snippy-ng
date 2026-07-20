from pathlib import Path
import sys

import click

from snippy_ng.utils.vcf_merge import (
    VcfMergeError,
    read_vcf_file_list,
    upset_svg_for_vcfs,
    upset_svg_from_merged_vcf_text,
)


@click.group("plot", context_settings={"show_default": True})
def plot() -> None:
    """VCF plotting utilities."""


@plot.command("upset", context_settings={"show_default": True})
@click.argument("vcfs", nargs=-1, type=click.Path(path_type=Path, exists=True, readable=True, dir_okay=False))
@click.option("--file-list", "file_lists", multiple=True, type=click.Path(path_type=Path, exists=True, readable=True, dir_okay=False), help="File with one VCF path per line; may be supplied more than once")
@click.option("--output", "-o", "output_svg", type=click.Path(path_type=Path, dir_okay=False), help="SVG output path (defaults to stdout)")
@click.option("--max-intersections", type=click.IntRange(min=1), default=40, help="Maximum intersections to show in the UpSet plot")
def upset(
    vcfs: tuple[Path, ...],
    file_lists: tuple[Path, ...],
    output_svg: Path | None,
    max_intersections: int,
) -> None:
    """Create an UpSet SVG from VCFs or merged VCF data on stdin."""
    input_vcfs = [*vcfs, *read_vcf_file_list(file_lists)]

    try:
        if input_vcfs:
            missing = [path for path in input_vcfs if not path.is_file()]
            if missing:
                raise click.UsageError(f"VCF does not exist: {missing[0]}")
            svg = upset_svg_for_vcfs(input_vcfs, max_intersections=max_intersections)
        else:
            stdin_text = sys.stdin.read()
            svg = upset_svg_from_merged_vcf_text(stdin_text, max_intersections=max_intersections)
    except VcfMergeError as exc:
        raise click.ClickException(str(exc)) from exc

    if output_svg is not None:
        output_svg.parent.mkdir(parents=True, exist_ok=True)
        output_svg.write_text(svg, encoding="utf-8")
        return
    click.echo(svg, nl=False)
