from pathlib import Path
from typing import Any
import click

from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={"show_default": True})
@add_snippy_global_options()
@click.argument(
    "config",
    required=True,
    type=click.File(mode="r")
)
@click.option(
    "--reference",
    "--ref",
    required=False,
    type=AbsolutePath(exists=True, readable=True),
)
@click.option("--cpus-per-sample", type=click.INT, default=1, help="Number of CPUs to allocate per sample")
@click.option("--core", type=click.FLOAT, default=0.95, help="Proportion of samples a site must be present in to be included in the core alignment (0.0-1.0)")
def multi(config: click.File, reference: Path | None, cpus_per_sample: int, core: float, outdir: Path, prefix: str, **context: Any):
    """
    Multi-sample SNP calling pipeline

    Example usage:

        $ snippy-ng multi samples.csv --ref reference.fasta
        $ snippy-ng gather | snippy-ng multi --ref reference.fasta - 

    """
    from snippy_ng import Context
    from snippy_ng.pipelines.common import load_or_prepare_reference
    from snippy_ng.pipelines import SnippyPipeline
    from snippy_ng.pipelines.multi import load_multi_config, run_multi_pipeline
    
    try:
        cfg = load_multi_config(config, reference)
    except Exception as e:
        raise click.ClickException(e)
    if reference is not None:
        cfg["reference"] = str(reference)
    # create reusable reference
    ref_stage = load_or_prepare_reference(
        reference_path=cfg["reference"],
        output_directory=Path(outdir) / 'reference',
    )
    ref_pipeline = SnippyPipeline(stages=[ref_stage])
    context["outdir"] = outdir
    run_ctx = Context(**context)
    ref_pipeline.run(run_ctx)

    snippy_reference_dir = ref_stage.output.reference.parent
    run_multi_pipeline(
        snippy_reference_dir=snippy_reference_dir,
        samples=cfg["samples"],
        outdir=outdir,
        prefix=prefix,
        tmpdir=context["tmpdir"],
        cpus=context["cpus"],
        ram=context["ram"],
        skip_check=context["skip_check"],
        check=context["check"],
        quiet=context["quiet"],
        create_missing=context["create_missing"],
        keep_incomplete=context["keep_incomplete"],
        cpus_per_sample=cpus_per_sample,
    )
    
    # core alignment
    from snippy_ng.pipelines.core import CorePipelineBuilder

    aln_pipeline = CorePipelineBuilder(
        snippy_dirs=[str(outdir / 'samples' / sample) for sample in cfg["samples"]],
        reference=snippy_reference_dir,
        core=core,
    ).build()
    core_outdir = Path(outdir) / 'core'
    core_outdir.mkdir(parents=True, exist_ok=True)
    context["outdir"] = core_outdir
    core_run_ctx = Context(**context)
    return aln_pipeline.run(core_run_ctx)



