import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Tuple
from snippy_ng.exceptions import PipelineExecutionError, SnippyError
from snippy_ng.logging import logger

def run_multi_pipeline(
    snippy_reference_dir: Path,
    samples: Dict[str, Any],
    config: Dict[str, Any],
) -> None:

    total_cpus = int(config["cpus"])
    cpus_per_sample = min(int(config["cpus_per_sample"]), total_cpus) if config["cpus_per_sample"] else total_cpus
    max_parallel = max(1, total_cpus // cpus_per_sample)


    # Minimal picklable config
    global_config = dict(config)
    global_config["cpus_per_sample"] = cpus_per_sample
    global_config['outdir'] = str(Path(config.get('outdir')).resolve())
    global_config['reference'] = str(snippy_reference_dir)

    jobs = [
        (sample_name, sample_cfg, global_config)
        for sample_name, sample_cfg in samples.items()
    ]

    ctx = mp.get_context("spawn")  # safest cross-platform
    failures = []

    with ProcessPoolExecutor(max_workers=max_parallel, mp_context=ctx) as executor:
        futures = {executor.submit(_run_one_sample, job): job[0] for job in jobs}

        for fut in as_completed(futures):
            sample = futures[fut]
            try:
                fut.result()
            except SnippyError as e:
                failures.append((sample, str(e)))
                logger.error(f"FAILED: {sample}: {e}")

    if failures:
        raise PipelineExecutionError("Some samples failed (check logs for details):\n" + "\n".join(f"- {s} -> {msg}" for s, msg in failures))
    

def _run_one_sample(job: Tuple[str, Dict[str, Any], Dict[str, Any]]) -> str:
    sample_name, sample_cfg, config = job

    # Import inside worker for clean spawn
    from snippy_ng.pipelines.asm import create_asm_pipeline
    from snippy_ng.pipelines.long import create_long_pipeline
    from snippy_ng.pipelines.short import create_short_pipeline
    import click
    
    sample_type = sample_cfg.get("type")
    del sample_cfg["type"]

    if sample_type == "short":
        reads = sample_cfg.get("reads")
        if not reads:
            # if reads not provided, expect left/right
            reads = [str(Path(r).resolve()) for r in (sample_cfg.get("left"), sample_cfg.get("right")) if r]
        pipeline = create_short_pipeline(
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
        pipeline = create_long_pipeline(
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
        pipeline = create_asm_pipeline(
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
    code = pipeline(
        skip_check=config["skip_check"],
        check=config["check"],
        outdir=outdir,
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )

    if code != 0:
            raise PipelineExecutionError(f"Sample '{sample_name}' failed with exit code {code}")

    return sample_name


def load_multi_config(config: Dict[str, Any]) -> Dict[str, Any]:
    config_path = config["config"]
    is_tsv = config_path.endswith(".tsv")
    is_csv = config_path.endswith(".csv")
    is_json = config_path.endswith(".json")

    if is_csv or is_tsv:
        if not config.get("reference"):
            raise ValueError("Reference must be provided option when using CSV/TSV config")
        import csv

        with open(config_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t" if is_tsv else ",")
            cfg = {"samples": {}}
            for row in reader:
                sample = row.get("sample")
                if not sample:
                    raise ValueError("Each row in the config file must have a 'sample' column with a unique sample name")
                if sample in cfg["samples"]:
                    raise ValueError(f"Duplicate sample name '{sample}' found in config file")
                row.pop("sample", None)
                cfg["samples"][sample] = {k: v for k, v in row.items() if v}
    elif is_json:
        import json

        with open(config_path) as f:
            cfg = json.load(f)
        if not config.get("reference") and "reference" not in cfg:
            raise ValueError("Reference must be provided in config file or as a command-line option")
        if "reference" in cfg and not config.get("reference"):
            config["reference"] = cfg["reference"]
    else:
        raise ValueError("Unsupported config file format. Use CSV, TSV, or JSON.")

    for name, row in cfg.get("samples", {}).items():
        if "type" not in row:
            raise ValueError(f"Sample '{name}' is missing required 'type' field")

    return cfg