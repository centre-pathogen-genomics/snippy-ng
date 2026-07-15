from pathlib import Path

import click

from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.utils.vcf import filter_variant_context_vcf


@click.command("context-filter", context_settings={"show_default": True})
@click.argument("input_vcf", type=AbsolutePath(exists=True, readable=True))
@click.option("--output", "-o", "output_vcf", required=True, type=click.Path(path_type=Path, dir_okay=False), help="Output VCF path")
@click.option("--max-local-snps", type=click.IntRange(min=0), default=0, help="Maximum SNPs allowed within the local SNP window; 0 disables")
@click.option("--local-snp-window", type=click.IntRange(min=0), default=0, help="Reference-base radius for the local SNP window; 0 disables")
@click.option("--min-snp-distance-to-indel", type=click.IntRange(min=0), default=0, help="Minimum SNP-to-indel distance; 0 disables")
@click.option("--min-snp-distance-to-breakpoint", type=click.IntRange(min=0), default=0, help="Minimum SNP-to-alignment-edge distance; 0 disables")
def context_filter(
    input_vcf: Path,
    output_vcf: Path,
    max_local_snps: int,
    local_snp_window: int,
    min_snp_distance_to_indel: int,
    min_snp_distance_to_breakpoint: int,
) -> None:
    """Apply local-context filters to a VCF."""
    if not any(
        [
            max_local_snps > 0 and local_snp_window > 0,
            min_snp_distance_to_indel > 0,
            min_snp_distance_to_breakpoint > 0,
        ]
    ):
        raise click.UsageError("At least one filter must be enabled")
    filter_variant_context_vcf(
        input_vcf,
        output_vcf,
        max_local_snps,
        local_snp_window,
        min_snp_distance_to_indel,
        min_snp_distance_to_breakpoint,
    )
