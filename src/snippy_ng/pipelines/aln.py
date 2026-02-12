from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.pipelines import SnippyPipeline, PipelineBuilder
from snippy_ng.stages.msa import CombineFastaFile, SoftCoreFilter
from typing import Optional
from pathlib import Path
from pydantic import Field


class AlnPipelineBuilder(PipelineBuilder):
    """Builder for alignment (MSA) pipeline."""
    snippy_dirs: list[Path] = Field(..., description="List of Snippy output directories")
    reference: Path = Field(..., description="Reference genome file")
    core: float = Field(default=0.95, description="Core genome threshold (0-1)")
    tmpdir: Optional[Path] = Field(default=None, description="Temporary directory")
    cpus: int = Field(default=1, description="Number of CPUs to use")
    ram: int = Field(default=8, description="RAM in GB")
    prefix: str = Field(default="core", description="Output file prefix")

    def build(self) -> SnippyPipeline:
        """Build and return the alignment pipeline."""
        stages = []

        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=self.reference
        )
        reference_file = setup.output.reference
        stages.append(setup)

        # Stage to combine FASTA files into a single alignment
        combine_stage = CombineFastaFile(
            snippy_dirs=self.snippy_dirs,
            reference=reference_file,
            tmpdir=self.tmpdir,
            cpus=self.cpus,
            ram=self.ram,
            prefix=self.prefix,
        )
        stages.append(combine_stage)

        # # Stage to filter the alignment to create core alignment
        filter_stage = SoftCoreFilter(
            aln=combine_stage.output.aln,
            core_threshold=self.core,
            tmpdir=self.tmpdir,
            cpus=self.cpus,
            ram=self.ram,
            prefix=self.prefix,
        )
        stages.append(filter_stage)

        return SnippyPipeline(stages=stages)