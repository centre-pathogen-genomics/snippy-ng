from pathlib import Path
from typing import List, Optional

import numpy as np

from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.dependencies import bedtools, bcftools, samtools

from pydantic import Field


def resolve_min_depth(min_depth: int, min_depth_file: Optional[Path]) -> int:
    """Return the effective minimum depth, preferring a dynamically computed value from file."""
    if min_depth_file is None:
        return min_depth
    return int(Path(min_depth_file).read_text().strip())


class MinDepthThresholdFromBamOutput(BaseOutput):
    """Output from dynamic minimum depth threshold computation."""

    depth_txt: TempPath = Field(..., description="Temporary per-base depth file")
    threshold_txt: Path = Field(
        ..., description="Computed minimum depth threshold (integer)"
    )


class MinDepthThresholdFromBam(BaseStage):
    """Compute a dynamic minimum depth threshold from per-base BAM depth.

    The threshold is max(min_depth, round(median(DEPTH) - 2 * MAD(DEPTH)))
    over covered positions (depth > 0), so unusually low-depth regions stand
    out from the bulk of the genome. If no positions are covered the median
    and MAD are treated as 0, so the threshold falls back to min_depth.
    """

    bam: Path = Field(..., description="Input BAM file")
    min_depth: int = Field(..., description="Minimum allowed threshold (lower bound)")

    _dependencies = [
        samtools,
    ]

    @property
    def output(self) -> MinDepthThresholdFromBamOutput:
        return MinDepthThresholdFromBamOutput(
            depth_txt=Path(f"{self.prefix}.mindepth.depth.txt"),
            threshold_txt=Path(f"{self.prefix}.mindepth.threshold.txt"),
        )

    def create_commands(self, ctx) -> List:
        return [
            self.shell_cmd(
                ["samtools", "depth", str(self.bam)],
                description="Generate per-base depth for dynamic minimum depth threshold",
                output_file=self.output.depth_txt,
            ),
            self.python_cmd(
                func=self.write_min_depth_threshold,
                args=[self.output.depth_txt, self.output.threshold_txt, self.min_depth],
                description="Compute minimum depth threshold from per-base depth (median - 2*MAD)",
            ),
        ]

    @staticmethod
    def write_min_depth_threshold(
        depth_txt: Path, threshold_txt: Path, min_depth: int
    ) -> None:
        """Write max(min_depth, round(median - 2*MAD)) over covered per-base depths."""
        import array

        # array('I') stores 4-byte unsigned ints; np.frombuffer views them as uint32
        depths = array.array("I")
        with open(depth_txt, "r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 3:
                    continue
                try:
                    depth = int(fields[2])
                except ValueError:
                    continue
                if depth > 0:
                    depths.append(depth)

        if not depths:
            median = 0.0
            mad = 0.0
        else:
            depth_array = np.frombuffer(depths, dtype=np.uint32)
            median = float(np.median(depth_array))
            mad = float(np.median(np.abs(depth_array - median)))
        threshold = max(min_depth, int(round(median - 2 * mad)))
        threshold_txt.write_text(f"{threshold}\n")
        print(
            f"Computed dynamic minimum depth threshold: {threshold}x "
            f"(median depth: {median:.0f}x, MAD: {mad:.0f}x)"
        )


class DepthBedsFromBamOutput(BaseOutput):
    """Output from combined depth BED generation."""

    zero_depth_bed: Path = Field(..., description="BED file with zero-depth regions")
    min_depth_bed: Path = Field(
        ..., description="BED file with regions to mask due to low depth"
    )


class DepthBedsFromBam(BaseStage):
    """Generate zero-depth and minimum-depth BED files from per-base depth."""

    bam: Path = Field(..., description="Input BAM file")
    min_depth: int = Field(..., description="Minimum depth threshold")
    min_depth_file: Optional[Path] = Field(
        default=None,
        description="File with a dynamically computed minimum depth threshold",
    )

    _dependencies = [
        bedtools,
        samtools,
    ]

    @property
    def output(self) -> DepthBedsFromBamOutput:
        return DepthBedsFromBamOutput(
            zero_depth_bed=Path(f"{self.prefix}.zerodepth.bed"),
            min_depth_bed=Path(f"{self.prefix}.mindepth.bed"),
        )

    def create_commands(self, ctx) -> List:
        min_depth = resolve_min_depth(self.min_depth, self.min_depth_file)
        depth_cmd = [
            "samtools",
            "depth",
            "-aa",
            str(self.bam),
        ]
        return [
            self.shell_pipe(
                [
                    self.shell_cmd(
                        depth_cmd,
                        description="Generate per-base depth for zero-depth mask",
                    ),
                    self.shell_cmd(
                        ["awk", '$3==0 {print $1"\\t"($2-1)"\\t"$2}'],
                        description="Extract zero-depth bases as BED intervals",
                    ),
                    self.shell_cmd(
                        ["bedtools", "merge", "-i", "-"],
                        description="Merge adjacent zero-depth BED intervals",
                    ),
                ],
                description="Generate zero-depth BED blocks",
                output_file=self.output.zero_depth_bed,
            ),
            self.shell_pipe(
                [
                    self.shell_cmd(
                        depth_cmd,
                        description="Generate per-base depth for minimum-depth mask",
                    ),
                    self.shell_cmd(
                        [
                            "awk",
                            f'$3>0 && $3<{min_depth} {{print $1"\\t"($2-1)"\\t"$2}}',
                        ],
                        description="Extract low-depth bases as BED intervals",
                    ),
                    self.shell_cmd(
                        ["bedtools", "merge", "-i", "-"],
                        description="Merge adjacent low-depth BED intervals",
                    ),
                ],
                description="Generate minimum-depth BED blocks",
                output_file=self.output.min_depth_bed,
            ),
        ]


class DepthMaskFromBedOutput(BaseOutput):
    """Output from applying a precomputed minimum-depth mask."""

    masked_fasta: Path = Field(
        ..., description="FASTA with low-depth regions (depth < min_depth) masked"
    )


class ApplyDepthMaskToFasta(BaseStage):
    """Apply a precomputed minimum-depth BED mask to a FASTA file."""

    fasta: Path = Field(..., description="Input FASTA file to be masked")
    mask_bed: Path = Field(
        ..., description="BED file with regions to mask due to low depth"
    )
    min_depth: int = Field(..., description="Minimum depth threshold")
    min_depth_file: Optional[Path] = Field(
        default=None,
        description="File with a dynamically computed minimum depth threshold",
    )
    mask_char: str = Field("N", description="Character to use for masking")

    _dependencies = [bedtools]

    @property
    def output(self) -> DepthMaskFromBedOutput:
        return DepthMaskFromBedOutput(
            masked_fasta=Path(f"{self.prefix}.mindepth_masked.fasta")
        )

    def create_commands(self, ctx) -> List:
        min_depth = resolve_min_depth(self.min_depth, self.min_depth_file)
        return [
            self.shell_cmd(
                [
                    "bedtools",
                    "maskfasta",
                    "-fi",
                    str(self.fasta),
                    "-bed",
                    str(self.mask_bed),
                    "-fo",
                    str(self.output.masked_fasta),
                    "-fullHeader",
                    "-mc",
                    self.mask_char,
                ],
                description=f"Apply min-depth mask (0 < depth < {min_depth})",
            )
        ]


class ApplyMaskOutput(BaseOutput):
    masked_fasta: Path = Field(
        ..., description="FASTA with user-supplied BED mask applied"
    )


class ApplyMask(BaseStage):
    """
    Masking stage that applies a supplied BED mask to a FASTA file.
    """

    fasta: Path = Field(..., description="Input FASTA file to be masked")
    mask_bed: Path = Field(..., description="BED file with regions to mask")
    mask_char: str = Field("X", description="Character to use for masking")

    _dependencies = [bedtools]

    @property
    def output(self) -> ApplyMaskOutput:
        return ApplyMaskOutput(masked_fasta=Path(f"{self.prefix}.masked.fasta"))

    def create_commands(self, ctx) -> List:
        """Apply mask to FASTA file using temporary copy"""
        temp_fasta = self.fasta.with_suffix(".tmp")

        return [
            # Copy input FASTA to temporary location
            self.shell_cmd(
                ["cp", str(self.fasta), str(temp_fasta)],
                description=f"Copy input FASTA to temporary location: {temp_fasta}",
            ),
            # Apply mask to temporary FASTA
            self.shell_cmd(
                [
                    "bedtools",
                    "maskfasta",
                    "-fi",
                    str(temp_fasta),
                    "-bed",
                    str(self.mask_bed),
                    "-fo",
                    str(self.output.masked_fasta),
                    "-fullHeader",
                    "-mc",
                    self.mask_char,
                ],
                description="Masking FASTA with provided BED file",
            ),
            # Clean up temporary file
            self.shell_cmd(
                ["rm", str(temp_fasta)], description="Remove temporary FASTA file"
            ),
        ]


class MaskMixedSitesOutput(BaseOutput):
    mixed_site_bed: TempPath = Field(
        ..., description="BED file of sites flagged with the MIXEDSITE VCF filter"
    )
    masked_fasta: Path = Field(..., description="FASTA with MIXEDSITE positions masked")


class MaskMixedSites(BaseStage):
    """Mask sites flagged with the MIXEDSITE VCF filter using bedtools maskfasta."""

    fasta: Path = Field(..., description="Consensus FASTA to mask")
    vcf: Path = Field(..., description="VCF containing MIXEDSITE-filtered records")
    mask_char: str = Field("n", description="Character to use for masking")

    _dependencies = [
        bcftools,
        bedtools,
    ]

    @property
    def output(self) -> MaskMixedSitesOutput:
        return MaskMixedSitesOutput(
            mixed_site_bed=Path(f"{self.prefix}.mixed_sites.bed"),
            masked_fasta=Path(f"{self.prefix}.mixed_masked.fasta"),
        )

    def create_commands(self, ctx) -> List:
        return [
            self.shell_pipe(
                [
                    self.shell_cmd(
                        [
                            "bcftools",
                            "query",
                            "-i",
                            'FILTER="MIXEDSITE"',
                            "-f",
                            r"%CHROM\t%POS0\t%REF\n",
                            str(self.vcf),
                        ],
                        description="Extract MIXEDSITE VCF records",
                    ),
                    self.shell_cmd(
                        ["awk", '{print $1"\\t"$2"\\t"($2+length($3))}'],
                        description="Convert MIXEDSITE records to BED intervals",
                    ),
                    self.shell_cmd(
                        ["sort", "-k1,1", "-k2,2n"],
                        description="Sort MIXEDSITE BED intervals",
                    ),
                    self.shell_cmd(
                        ["bedtools", "merge", "-i", "-"],
                        description="Merge adjacent MIXEDSITE BED intervals",
                    ),
                ],
                description="Create BED file for MIXEDSITE positions",
                output_file=self.output.mixed_site_bed,
            ),
            self.shell_cmd(
                [
                    "bedtools",
                    "maskfasta",
                    "-fi",
                    str(self.fasta),
                    "-bed",
                    str(self.output.mixed_site_bed),
                    "-fo",
                    str(self.output.masked_fasta),
                    "-fullHeader",
                    "-mc",
                    self.mask_char,
                ],
                description="Mask MIXEDSITE positions in consensus FASTA",
            ),
        ]
