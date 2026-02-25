from pathlib import Path
from typing import Optional, Union
from pydantic import Field
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.reporting import EpiReport, Context


class ReportPipelineBuilder(PipelineBuilder):
    """Builder for reporting pipeline."""
    tree: Path = Field(..., description="Phylogenetic tree file in Newick format")
    preprocess_tree: bool = Field(default=True, description="Whether to pre-process the NEWICK tree by rooting at midpoint and ladderizing before rendering the report")
    title: str = Field(default="Snippy-NG Report", description="Title for the report")
    metadata: Optional[Union[Path, str]] = Field(default=None, description="Metadata file in JSON format or a JSON string")
    logs: Optional[Path] = Field(default=None, description="Log file to include in the report")
    

    def build(self) -> SnippyPipeline:
        """Build and return the tree building pipeline."""
        stages = []
        context: Context = {
                "NEWICK": self.tree, 
                "REPORT_NAME": self.title, 
                "METADATA_JSON": self.metadata if self.metadata else None,
                "LOGS": self.logs if self.logs else None,
            }
        report_stage = EpiReport(
            context=context,
            preprocess_tree=self.preprocess_tree,
        )
        stages.append(report_stage)

        return SnippyPipeline(stages=stages)