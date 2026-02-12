from pathlib import Path
from typing import Optional
from pydantic import Field
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.trees import IQTreeBuildTree


class TreePipelineBuilder(PipelineBuilder):
    """Builder for phylogenetic tree building pipeline."""
    aln: str = Field(..., description="Multiple sequence alignment file")
    model: str = Field(default="GTR+G", description="Substitution model for IQ-TREE")
    bootstrap: int = Field(default=1000, description="Number of bootstrap replicates")
    fconst: Optional[str] = Field(default=None, description="Frequency constants for ascertainment bias correction")
    tmpdir: Optional[Path] = Field(default=None, description="Temporary directory")
    cpus: int = Field(default=1, description="Number of CPUs to use")
    ram: int = Field(default=8, description="RAM in GB")
    prefix: str = Field(default="tree", description="Output file prefix")

    def build(self) -> SnippyPipeline:
        """Build and return the tree building pipeline."""
        stages = []

        # Stage to build a phylogenetic tree using IQ-TREE
        iqtree_stage = IQTreeBuildTree(
            aln=self.aln,
            model=self.model,
            bootstrap=self.bootstrap,
            fconst=self.fconst,
            tmpdir=self.tmpdir,
            cpus=self.cpus,
            ram=self.ram,
            prefix=self.prefix,
        )
            
        stages.append(iqtree_stage)

        return SnippyPipeline(stages=stages)