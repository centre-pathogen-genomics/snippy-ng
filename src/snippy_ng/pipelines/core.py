from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.pipelines import SnippyPipeline, PipelineBuilder
from snippy_ng.stages.core import CombineFastaFile, SoftCoreFilter
from snippy_ng.stages.stats import AlignmentAlignedPercentage
from pathlib import Path
from pydantic import Field


class CorePipelineBuilder(PipelineBuilder):
    """Builder for alignment (MSA) pipeline."""
    snippy_dirs: list[Path] = Field(..., description="List of Snippy output directories")
    reference: Path = Field(..., description="Reference genome file")
    core: float = Field(default=0.95, description="Core genome threshold (0-1)")
    prefix: str = Field(default="core", description="Output file prefix")

    def build(self) -> SnippyPipeline:
        """Build and return the alignment pipeline."""
        stages = []

        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=self.reference,
            output_directory=Path("reference"),
        )
        reference_file = setup.output.reference
        reference_id = ReferenceMetadata(setup.output.metadata).prefix or "reference"
        stages.append(setup)

        # Stage to combine FASTA files into a single alignment
        combine_stage = CombineFastaFile(
            snippy_dirs=self.snippy_dirs,
            reference=reference_file,
            reference_id=reference_id,
            prefix=self.prefix,
        )
        stages.append(combine_stage)

        alignment_stats = AlignmentAlignedPercentage(
            alignment=combine_stage.output.aln,
            prefix=self.prefix,
        )
        stages.append(alignment_stats)

        # Stage to filter the alignment to create core alignment
        filter_stage = SoftCoreFilter(
            aln=combine_stage.output.aln,
            core_threshold=self.core,
            prefix=self.prefix,
        )
        stages.append(filter_stage)

        return SnippyPipeline(stages=stages)
