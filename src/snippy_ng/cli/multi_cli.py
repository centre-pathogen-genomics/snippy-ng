from pathlib import Path
from typing import Dict, Any, Tuple
import click


from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(cls=CommandWithGlobals, context_settings={"show_default": True})
@add_snippy_global_options()
@click.option(
    "--config",
    required=True,
    type=click.Path(exists=True, resolve_path=True, readable=True),
)
@click.option("--cpus-per-sample", type=click.INT, default=1, show_default=True)
@click.option(
    "--reference",
    "--ref",
    required=True,
    type=click.Path(exists=True, resolve_path=True, readable=True),
)
@click.option("--core", type=click.FLOAT, default=0.95, help="Proportion of samples a site must be present in to be included in the core alignment (0.0-1.0)")
def multi(**config):
    from snippy_ng.cli.utils.pipeline_runner import run_snippy_pipeline
    from snippy_ng.pipelines.common import prepare_reference
    import yaml
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    # Load YAML
    with open(config["config"], "r") as f:
        try:
            cfg = yaml.safe_load(f) or {}
        except yaml.YAMLError as e:
            raise click.UsageError(f"Error parsing configuration file: {e}")
        
    # create reusable reference
    reference_outdir = config["outdir"] / "reference"
    code = run_snippy_pipeline(
        stages=[
            prepare_reference(
                reference_path=config["reference"],
                output_directory=reference_outdir
            )
        ],
        skip_check=config["skip_check"],
        check=config["check"],
        outdir=Path("."),
        quiet=config["quiet"],
        continue_last_run=config["continue_last_run"],
        keep_incomplete=config["keep_incomplete"],
    )
    if code != 0:
        raise click.ClickException("Reference preparation failed, aborting multi-sample run.")

    total_cpus = int(config["cpus"])
    cpus_per_sample = int(config["cpus_per_sample"])
    max_parallel = max(1, total_cpus // cpus_per_sample)


    # Minimal picklable config
    global_config = dict(config)
    global_config["cpus_per_sample"] = cpus_per_sample
    global_config['outdir'] = str(Path(config.get('outdir')).resolve())
    global_config['reference'] = str(reference_outdir.resolve())

    jobs = [
        (sample_name, sample_cfg, global_config)
        for sample_name, sample_cfg in cfg['samples'].items()
    ]

    ctx = mp.get_context("spawn")  # safest cross-platform
    failures = []

    with ProcessPoolExecutor(max_workers=max_parallel, mp_context=ctx) as executor:
        futures = {executor.submit(_run_one_sample, job): job[0] for job in jobs}

        for fut in as_completed(futures):
            sample = futures[fut]
            try:
                fut.result()
            except Exception as e:
                failures.append((sample, str(e)))
                click.echo(f"FAILED: {sample}: {e}", err=True)

    if failures:
        raise click.ClickException(
            "Some samples failed:\n" + "\n".join(f"- {s}: {msg}" for s, msg in failures)
        )
    
    # core alignment
    from snippy_ng.pipelines.aln import create_aln_pipeline_stages

    stages = create_aln_pipeline_stages(
        snippy_dirs=[str((Path(config["outdir"]) / 'samples' / sample).resolve()) for sample, *_ in jobs],
        reference=global_config["reference"],
        core=config["core"],
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    )
    outdir = Path(config['outdir']) / 'core'
    outdir.mkdir(parents=True, exist_ok=True)
    return run_snippy_pipeline(
        stages,
        skip_check=config['skip_check'],
        check=config['check'],
        outdir=outdir,
        quiet=config['quiet'],
        continue_last_run=config['continue_last_run'],
        keep_incomplete=config['keep_incomplete'],
    )



def _run_one_sample(job: Tuple[str, Dict[str, Any], Dict[str, Any]]) -> str:
    sample_name, sample_cfg, config = job

    # Import inside worker for clean spawn
    from snippy_ng.pipelines.asm import create_asm_pipeline_stages
    from snippy_ng.pipelines.long import create_long_pipeline_stages
    from snippy_ng.pipelines.short import create_short_pipeline_stages
    from snippy_ng.cli.utils.pipeline_runner import run_snippy_pipeline
    import click
    
    sample_type = sample_cfg.get("type")
    del sample_cfg["type"]

    if sample_type == "short":
        reads = sample_cfg.get("reads")
        if not reads:
            # if reads not provided, expect left/right
            reads = [str(Path(r).resolve()) for r in (sample_cfg.get("left"), sample_cfg.get("right")) if r]
        stages = create_short_pipeline_stages(
            reference=config["reference"],
            reads=reads,
            bam=str(Path(sample_cfg.get("bam")).resolve()) if sample_cfg.get("bam") else None,
            prefix=config["prefix"],
            tmpdir=config["tmpdir"],
            cpus=config["cpus_per_sample"],
            ram=config["ram"],
            **{k: v for k, v in sample_cfg.items() if k not in ["left", "right", "bam", "reads"]},
        )

    elif sample_type == "long":
        stages = create_long_pipeline_stages(
            reference=config["reference"],
            reads=str(Path(sample_cfg.get("reads")).resolve()) if sample_cfg.get("reads") else None,
            bam=str(Path(sample_cfg.get("bam")).resolve()) if sample_cfg.get("bam") else None,
            prefix=config["prefix"],
            tmpdir=config["tmpdir"],
            cpus=config["cpus_per_sample"],
            ram=config["ram"],
            **{k: v for k, v in sample_cfg.items() if k not in ["reads", "bam"]},
        )

    elif sample_type == "asm":
        sample_cfg['assembly'] = str(Path(sample_cfg.get("assembly")).resolve())
        stages = create_asm_pipeline_stages(
            reference=config["reference"],
            prefix=config["prefix"],
            tmpdir=config["tmpdir"],
            cpus=config["cpus_per_sample"],
            ram=config["ram"],
            **sample_cfg,
        )

    else:
        raise click.UsageError(
            f"Unknown sample type '{sample_type}' for sample '{sample_name}'"
        )

    # Per-sample output directory to prevent collisions
    outdir = (Path(config["outdir"]) / 'samples' / sample_name).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    
    # run_snippy_pipeline sets the working dir to outdir
    code = run_snippy_pipeline(
        stages=stages,
        skip_check=config["skip_check"],
        check=config["check"],
        outdir=outdir,
        quiet=config["quiet"],
        continue_last_run=config["continue_last_run"],
        keep_incomplete=config["keep_incomplete"],
    )

    if code != 0:
        raise RuntimeError(f"Sample '{sample_name}' failed with exit code {code}")

    return sample_name
