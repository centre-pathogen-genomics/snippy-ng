from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.pipelines import SnippyPipeline, PipelineBuilder
from snippy_ng.stages.alignment_filter import CheckAlignmentClustersByPipelineType, FilterAlignmentByAlignedPercentage
from snippy_ng.stages.core import CombineFastaFile, DistleDistanceMatrix, SoftCoreFilter
from snippy_ng.stages.stats import AlignmentAlignedPercentage
from pathlib import Path
from pydantic import Field


class CorePipelineBuilder(PipelineBuilder):
    """Builder for alignment (MSA) pipeline."""
    snippy_dirs: list[Path] = Field(..., description="List of Snippy output directories")
    reference: Path = Field(..., description="Reference genome file")
    core: float = Field(default=0.95, description="Core genome threshold (0-1)")
    inclusion_threshold: float = Field(default=0.20, ge=0.0, le=1.0, description="Posterior probability threshold for retaining membership in the retained alignment percentage clusters")
    snp_distance_format: str = Field(default="tabular", description="Output format for pairwise snip distance matrix")
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
            format=self.snp_distance_format,
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
        outputs_to_keep.extend(alignment_filter.output.paths)

        technical_cluster_check = CheckAlignmentClustersByPipelineType(
            filter_stats=alignment_filter.output.filter_stats,
            qc_files=[
                qc_file
                for sample_dir in self.snippy_dirs
                for qc_file in sorted(sample_dir.glob("*.qc.tsv"))
            ],
            prefix=self.prefix,
        )
        stages.append(technical_cluster_check)
        outputs_to_keep.extend(technical_cluster_check.output.paths)

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
            format=self.snp_distance_format,
        )
        stages.append(soft_core_distances)
        outputs_to_keep.extend(soft_core_distances.output.paths)

        return SnippyPipeline(stages=stages, outputs_to_keep=outputs_to_keep)
