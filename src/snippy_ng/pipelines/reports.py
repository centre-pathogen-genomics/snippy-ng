from pathlib import Path
from typing import Optional, Union
from pydantic import Field
import snippy_ng.pipelines as pipelines
from snippy_ng.pipelines import PipelineBuilder
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.stages.downsample import SamtoolsDownsampleAlignment
from snippy_ng.stages.reporting import SampleReport, TreeReport, Context


class TreeReportPipelineBuilder(PipelineBuilder):
    """Builder for reporting pipeline."""
    tree: Path = Field(..., description="Phylogenetic tree file in Newick format")
    mid_point_root: bool = Field(default=False, description="Mid-point root the tree in the report")
    ladderize: bool = Field(default=False, description="Ladderize the tree in the report")
    title: str = Field(default="Snippy-NG Report", description="Title for the report")
    metadata: Optional[Union[Path,str]] = Field(default=None, description="Metadata JSON string or file path to include in the report")
    color_by_column: Optional[str] = Field(default=None, description="Column name in the metadata to color the tree by")
    logs: Optional[Path] = Field(default=None, description="Log file to include in the report")
    

    def build(self):
        """Build and return the tree building pipeline."""
        stages = []
        context: Context = {
                "NEWICK": self.tree, 
                "REPORT_NAME": self.title, 
                "METADATA": lambda: self.metadata if self.metadata else None,
                "COLOR_BY_COLUMN": self.color_by_column if self.color_by_column else "",
                "LOGS": self.logs if self.logs else None,
            }
        report_stage = TreeReport(
            context=context,
            prefix=self.prefix,
            mid_point_root=self.mid_point_root,
            ladderize=self.ladderize,
        )
        stages.append(report_stage)

        return pipelines.SnippyPipeline(stages=stages)


class SampleReportPipelineBuilder(PipelineBuilder):
    """Builder for per-sample variant and alignment reports."""

    vcf: Path = Field(..., description="VCF file used to build the variant table")
    alignment: Optional[Path] = Field(default=None, description="Optional BAM/CRAM alignment to window around variants")
    reference: Optional[Path] = Field(default=None, description="Optional reference genome (FASTA or GenBank) or prepared reference directory")
    title: str = Field(default="Snippy-NG Sample Report", description="Title for the report")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override")
    variant_scope: str = Field(default="all", description="Variant scope to include: pass or all")
    window_size: int = Field(default=100, description="Base pairs of context to embed around each variant")
    downsample: Optional[float] = Field(default=None, gt=0, lt=1, description="Optional fraction of alignments to keep before report windowing")

    def build(self):
        stages=[]
        
        if self.alignment and not self.reference:
            raise ValueError("--reference is required when --alignment is provided")
        if self.downsample is not None and not self.alignment:
            raise ValueError("alignment is required when downsample is provided")

        setup = None
        if self.alignment and self.reference:
            setup = load_or_prepare_reference(
                reference_path=self.reference,
                output_directory=Path("reference"),
            )
            stages.append(setup)

        report_alignment = self.alignment
        if self.alignment and self.downsample is not None:
            downsample_stage = SamtoolsDownsampleAlignment(
                alignment=self.alignment,
                fraction=self.downsample,
                reference=setup.output.reference if setup else None,
                prefix=self.prefix,
            )
            stages.append(downsample_stage)
            report_alignment = downsample_stage.output.bam

        report_stage = SampleReport(
            vcf=self.vcf,
            alignment=report_alignment,
            reference=setup.output.reference if setup else None,
            reference_index=setup.output.reference_index if setup else None,
            title=self.title,
            sample_name=self.sample_name,
            variant_scope=self.variant_scope,
            window_size=self.window_size,
            prefix=self.prefix,
        )
        stages.append(report_stage)

        return pipelines.SnippyPipeline(
            stages=stages,
            outputs_to_keep=[report_stage.output.rendered],
        )
