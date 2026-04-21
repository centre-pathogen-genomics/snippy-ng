from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.pipelines import SnippyPipeline, PipelineBuilder
from snippy_ng.stages.core import CombineFastaFile, DistleDistanceMatrix, FilterAlignmentByAlignedPercentage, SoftCoreFilter
from snippy_ng.stages.stats import AlignmentAlignedPercentage
from pathlib import Path
from pydantic import Field


class CorePipelineBuilder(PipelineBuilder):
    """Builder for alignment (MSA) pipeline."""
    snippy_dirs: list[Path] = Field(..., description="List of Snippy output directories")
    reference: Path = Field(..., description="Reference genome file")
    core: float = Field(default=0.95, description="Core genome threshold (0-1)")
    inclusion_threshold: float = Field(default=0.1, description="Posterior probability threshold for retaining membership in the main alignment percentage cluster")
    prefix: str = Field(default="core", description="Output file prefix")

    def build(self) -> SnippyPipeline:
        """Build and return the alignment pipeline."""
        stages = []
        outputs_to_keep = []

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
        outputs_to_keep.extend(combine_stage.output.paths)

        full_distances = DistleDistanceMatrix(
            aln=combine_stage.output.aln,
        )
        stages.append(full_distances)
        outputs_to_keep.extend(full_distances.output.paths)

        alignment_stats = AlignmentAlignedPercentage(
            alignment=combine_stage.output.aln,
            prefix=self.prefix,
        )
        stages.append(alignment_stats)
        outputs_to_keep.extend(alignment_stats.output.paths)

        alignment_filter = FilterAlignmentByAlignedPercentage(
            aln=combine_stage.output.aln,
            alignment_stats=alignment_stats.output.aligned_tsv,
            prefix=self.prefix,
            inclusion_threshold=self.inclusion_threshold,
        )
        stages.append(alignment_filter)

        # Stage to filter the alignment to create core alignment
        filter_stage = SoftCoreFilter(
            aln=alignment_filter.output.filtered_aln,
            core_threshold=self.core,
            prefix=self.prefix,
        )
        stages.append(filter_stage)
        outputs_to_keep.extend(filter_stage.output.paths)

        soft_core_distances = DistleDistanceMatrix(
            aln=filter_stage.output.soft_core,
        )
        stages.append(soft_core_distances)
        outputs_to_keep.extend(soft_core_distances.output.paths)

        return SnippyPipeline(stages=stages, outputs_to_keep=outputs_to_keep)
