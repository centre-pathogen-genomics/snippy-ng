from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.pipelines import SnippyPipeline
from snippy_ng.stages.msa import CombineFastaFile, SoftCoreFilter
from typing import Optional
from pathlib import Path

def create_aln_pipeline(
    snippy_dirs: list[Path],
    reference: Path,
    core: float,
    tmpdir: Optional[str],
    cpus: int,
    ram: int,
    prefix: str = "core"
) -> SnippyPipeline:

    stages = []

    # Setup reference (load existing or prepare new)
    setup = load_or_prepare_reference(
        reference_path=reference
    )
    reference_file = setup.output.reference
    stages.append(setup)

    # Stage to combine FASTA files into a single alignment
    combine_stage = CombineFastaFile(
        snippy_dirs=snippy_dirs,
        reference=reference_file,
        tmpdir=tmpdir,
        cpus=cpus,
        ram=ram,
        prefix=prefix,
    )
    stages.append(combine_stage)

    # # Stage to filter the alignment to create core alignment
    filter_stage = SoftCoreFilter(
        aln=combine_stage.output.aln,
        core_threshold=core,
        tmpdir=tmpdir,
        cpus=cpus,
        ram=ram,
        prefix=prefix,
    )
    stages.append(filter_stage)

    return SnippyPipeline(stages=stages)