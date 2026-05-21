from pathlib import Path
from typing import Optional, Union
from pydantic import Field
import snippy_ng.pipelines as pipelines
from snippy_ng.pipelines import PipelineBuilder
from snippy_ng.pipelines.common import load_or_prepare_reference
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
    alignment: Path = Field(..., description="BAM/CRAM alignment to window around variants")
    reference: Path = Field(..., description="Reference genome (FASTA or GenBank) or prepared reference directory")
    title: str = Field(default="Snippy-NG Sample Report", description="Title for the report")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override")
    variant_scope: str = Field(default="pass", description="Variant scope to include: pass or all")
    window_size: int = Field(default=100, description="Base pairs of context to embed around each variant")

    def build(self):
        stages=[]
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=self.reference,
            output_directory=Path("reference"),
        )
        stages.append(setup)

        report_stage = SampleReport(
            vcf=self.vcf,
            alignment=self.alignment,
            reference=setup.output.reference,
            reference_index=setup.output.reference_index,
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
