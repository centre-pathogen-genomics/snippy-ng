from pathlib import Path
from typing import Optional, Union
from pydantic import Field
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.reporting import EpiReport, Context


class ReportPipelineBuilder(PipelineBuilder):
    """Builder for reporting pipeline."""
    tree: Path = Field(..., description="Phylogenetic tree file in Newick format")
    mid_point_root: bool = Field(default=False, description="Mid-point root the tree in the report")
    ladderize: bool = Field(default=False, description="Ladderize the tree in the report")
    title: str = Field(default="Snippy-NG Report", description="Title for the report")
    metadata: Optional[Union[Path,str]] = Field(default=None, description="Metadata JSON string or file path to include in the report")
    color_by_column: Optional[str] = Field(default=None, description="Column name in the metadata to color the tree by")
    logs: Optional[Path] = Field(default=None, description="Log file to include in the report")
    

    def build(self) -> SnippyPipeline:
        """Build and return the tree building pipeline."""
        stages = []
        context: Context = {
                "NEWICK": self.tree, 
                "REPORT_NAME": self.title, 
                "METADATA": lambda: self.metadata if self.metadata else None,
                "COLOR_BY_COLUMN": self.color_by_column if self.color_by_column else "",
                "LOGS": self.logs if self.logs else None,
            }
        report_stage = EpiReport(
            context=context,
            prefix=self.prefix,
            mid_point_root=self.mid_point_root,
            ladderize=self.ladderize,
        )
        stages.append(report_stage)

        return SnippyPipeline(stages=stages)