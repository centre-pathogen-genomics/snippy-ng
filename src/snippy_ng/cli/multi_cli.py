from pathlib import Path
import click

from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={"show_default": True})
@add_snippy_global_options()
@click.argument(
    "config",
    required=True,
    type=click.Path(exists=True, resolve_path=True, readable=True),
)
@click.option(
    "--reference",
    "--ref",
    required=False,
    type=click.Path(exists=True, resolve_path=True, readable=True),
)
@click.option("--cpus-per-sample", type=click.INT, help="Number of CPUs to allocate per sample")
@click.option("--core", type=click.FLOAT, default=0.95, help="Proportion of samples a site must be present in to be included in the core alignment (0.0-1.0)")
def multi(**config):
    """
    Multi-sample SNP calling pipeline

    Example usage:

        snippy-ng multi samples.csv --ref reference.fasta
    """
    from snippy_ng.pipelines.common import load_or_prepare_reference
    from snippy_ng.pipelines import SnippyPipeline
    from snippy_ng.pipelines.multi import load_multi_config, run_multi_pipeline
    
    try:
        cfg = load_multi_config(config)
    except Exception as e:
        raise click.ClickException(e)
    if config["reference"] is not None:
        cfg["reference"] = str(config["reference"])
    # create reusable reference
    ref_stage = load_or_prepare_reference(
        reference_path=cfg["reference"],
        output_directory=Path(config["outdir"]) / 'reference',
    )
    ref_pipeline = SnippyPipeline(stages=[ref_stage])
    ref_pipeline(
        skip_check=config["skip_check"],
        check=config["check"],
        cwd=config["outdir"],
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )

    snippy_reference_dir = ref_stage.output.reference.parent
    run_multi_pipeline(
        snippy_reference_dir=snippy_reference_dir,
        samples=cfg["samples"],
        config=config,
    )
    
    # core alignment
    from snippy_ng.pipelines.aln import create_aln_pipeline

    aln_pipeline = create_aln_pipeline(
        snippy_dirs=[str((Path(config["outdir"]) / 'samples' / sample).resolve()) for sample in cfg["samples"]],
        reference=snippy_reference_dir,
        core=config["core"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    )
    outdir = Path(config['outdir']) / 'core'
    outdir.mkdir(parents=True, exist_ok=True)
    aln_pipeline(
        skip_check=config['skip_check'],
        check=config['check'],
        cwd=outdir,
        quiet=config['quiet'],
        create_missing=config['create_missing'],
        keep_incomplete=config['keep_incomplete'],
    )



