import json
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Optional
import os

import gradio as gr

from snippy_ng.cli.multi_cli import run_multi_config
from snippy_ng.utils.gather import gather_samples_config


def _as_path(value: Any) -> Optional[Path]:
    if value is None:
        return None
    if isinstance(value, dict):
        path = value.get("path") or value.get("name")
        return Path(path) if path else None
    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        return Path(value)
    path = getattr(value, "path", None) or getattr(value, "name", None)
    return Path(path) if path else None


def _paths_from_files(values: Any) -> list[Path]:
    if values is None:
        return []
    if not isinstance(values, list):
        values = [values]
    paths: list[Path] = []
    for value in values:
        path = _as_path(value)
        if path is not None:
            paths.append(path)
    return paths


def _original_filename(value: Any) -> Optional[str]:
    if isinstance(value, dict):
        name = value.get("orig_name") or value.get("name") or value.get("path")
    else:
        name = getattr(value, "orig_name", None) or getattr(value, "name", None) or getattr(value, "path", None)
    return Path(name).name if name else None


def _unique_destination(directory: Path, filename: str) -> Path:
    destination = directory / filename
    if not destination.exists():
        return destination
    stem = destination.stem
    suffix = destination.suffix
    for index in range(1, 1000):
        candidate = directory / f"{stem}-{index}{suffix}"
        if not candidate.exists():
            return candidate
    raise gr.Error(f"Could not create a unique upload staging path for {filename}.")


def _stage_uploaded_file(value: Any, directory: Path) -> Optional[Path]:
    source = _as_path(value)
    if source is None:
        return None
    filename = _original_filename(value) or source.name
    directory.mkdir(parents=True, exist_ok=True)
    destination = _unique_destination(directory, filename)
    try:
        os.symlink(source, destination)
    except OSError:
        shutil.copy2(source, destination)
    return destination


def _stage_uploaded_files(values: Any, directory: Path) -> list[Path]:
    if values is None:
        return []
    if not isinstance(values, list):
        values = [values]
    paths: list[Path] = []
    for value in values:
        path = _stage_uploaded_file(value, directory)
        if path is not None:
            paths.append(path)
    return paths


def _path_from_text(value: str | None) -> Optional[Path]:
    if not value or not value.strip():
        return None
    return Path(value.strip()).expanduser().absolute()


def _paths_from_text(value: str | None) -> list[Path]:
    if not value:
        return []
    return [Path(line.strip()).expanduser().absolute() for line in value.splitlines() if line.strip()]


def _gradio_temp_base() -> Path:
    return Path(os.environ.get("GRADIO_TEMP_DIR") or tempfile.gettempdir()).absolute()


def _validated_outdir_name(outdir_name: str | None) -> str:
    name = (outdir_name or _default_outdir_name()).strip()
    path = Path(name)
    if (
        "/" in name
        or "\\" in name
        or any(char.isspace() for char in name)
        or path.is_absolute()
        or len(path.parts) != 1
        or path.name in {"", ".", ".."}
    ):
        raise gr.Error("Output directory name must be a single directory name without spaces or path separators.")
    return path.name


def _outdir_from_name(outdir_name: str | None, temp_output: bool) -> Path:
    name = _validated_outdir_name(outdir_name)
    if temp_output:
        return Path(
            tempfile.mkdtemp(
                prefix=f"{name}-",
                dir=_gradio_temp_base(),
            )
        ).absolute()
    return (Path.cwd() / name).absolute()


def _default_outdir_name() -> str:
    return f"snippy-ng-{datetime.now():%Y%m%d-%H%M%S}"


def select_download_file(filename: str | None):
    if not filename:
        return gr.update(value=None, visible=False)
    path = Path(filename)
    if not path.is_file():
        return gr.update(value=None, visible=False)
    return gr.update(value=str(path), visible=True)


def _cleanup_temp_tarball(tarball: Any):
    path = _as_path(tarball)
    if path is not None and path.is_file() and path.parent.name.startswith("snippy-ng-tarball-"):
        shutil.rmtree(path.parent, ignore_errors=True)
    return gr.update(value=None, visible=False)


def create_output_tarball(output_dir: str | Path | None, previous_tarball: Any = None):
    _cleanup_temp_tarball(previous_tarball)
    if not output_dir:
        return gr.update(value=None, visible=False)
    path = Path(output_dir)
    if not path.is_dir():
        return gr.update(value=None, visible=False)
    archive_dir = Path(tempfile.mkdtemp(prefix="snippy-ng-tarball-", dir=_gradio_temp_base()))
    tarball = shutil.make_archive(
        str(archive_dir / path.name),
        "gztar",
        root_dir=path.parent,
        base_dir=path.name,
    )
    return gr.update(value=tarball, visible=True)


def _disable_run_button():
    return gr.update(interactive=False)


def _enable_run_button():
    return gr.update(interactive=True)


def run_snippy_multi_gui(
    reference_file: Any,
    reference_path: str,
    sample_files: Any,
    sample_paths: str,
    outdir_name: str,
    temp_output: bool,
    prefix: str,
    cpus: int,
    cpus_per_sample: int,
    ram: int,
    max_depth: int,
    aggressive_ids: bool,
    exclude_name_regex: str,
    core: float,
    inclusion_threshold: float,
    stop_on_failure: bool,
    skip_check: bool,
    force: bool,
) -> tuple[Any, Any]:
    prefix = (prefix or "snippy").strip() or "snippy"

    try:
        outdir_path = _outdir_from_name(outdir_name, bool(temp_output))
        if outdir_path.exists() and any(outdir_path.iterdir()) and not force:
            raise ValueError(f"Output directory '{outdir_path}' already exists and is not empty. Enable force overwrite to continue.")

        uploaded_reference = _stage_uploaded_file(reference_file, outdir_path / ".gui-uploads" / "reference")
        reference = uploaded_reference or _path_from_text(reference_path)
        if reference is None:
            raise gr.Error("Upload a reference FASTA or GenBank file, or enter a server reference path.")

        inputs = [
            *_stage_uploaded_files(sample_files, outdir_path / ".gui-uploads" / "samples"),
            *_paths_from_text(sample_paths),
        ]
        if not inputs:
            raise gr.Error("Upload sample files or enter one or more server sample paths.")

        config = gather_samples_config(
            inputs,
            max_depth=max(0, int(max_depth)),
            aggressive_ids=bool(aggressive_ids),
            exclude_name_regex=exclude_name_regex.strip() or None,
            reference=reference,
        )
        reference = config.get("reference")
        if not reference:
            raise ValueError("Reference genome is required.")
        samples = config.get("samples") or {}
        if not samples:
            raise ValueError("No samples were found. Check the selected files, paths, depth, and exclude pattern.")

        result = run_multi_config(
            config,
            reference=Path(reference),
            cpus_per_sample=max(1, int(cpus_per_sample)),
            core=float(core),
            inclusion_threshold=float(inclusion_threshold),
            stop_on_failure=bool(stop_on_failure),
            outdir=outdir_path,
            prefix=prefix,
            context={
                "outdir": outdir_path,
                "log_path": Path("LOG.txt"),
                "cpus": max(1, int(cpus)),
                "ram": max(1, int(ram)),
                "skip_check": bool(skip_check),
                "force": bool(force),
            },
        )
        successful_samples = result["successful_samples"]
        failed_samples = result["failed_samples"]
    except Exception as exc:
        status = f"Error: {str(exc)}"
        return gr.update(value=status, visible=True), None

    status = [
            "Snippy-NG multi run completed.",
            f"Output directory: {outdir_path}",
            f"Samples completed: {len(successful_samples)}",
            f"Samples failed: {len(failed_samples)}",
        ]
    if failed_samples:
        status.append(f"Error details for failed samples: {json.dumps(failed_samples, indent=2)}")
    return gr.update(value="\n".join(status), visible=True), str(outdir_path)


def create_app(temp_output: bool = False, server_paths: bool = True, max_cpus: int | None = None) -> gr.Blocks:
    """
    Create and return the Gradio interface for the Snippy-NG GUI.
    """
    detected_cpus = os.cpu_count() or 1
    max_cpus = max(1, int(max_cpus)) if max_cpus is not None else detected_cpus
    default_cpus = min(detected_cpus, max_cpus)
    default_cpus_per_sample = min(8, default_cpus)

    with gr.Blocks(title="Snippy-NG GUI") as app:
        input_message = (
            "Upload files from your computer or enter paths available on the server. "
            if server_paths
            else "Upload files from your computer. "
        )
        gr.Markdown(
            f"# Snippy-NG GUI\n"
            f"{input_message}"
            "Snippy-NG will create a multi-sample config, run `snippy-ng multi`, "
            "and make the output files available for download. "
            "[Read the docs](https://cpg.org.au/snippy-ng/)."
        )

        if server_paths:
            with gr.Tabs():
                with gr.Tab("Upload files"):
                    reference_file = gr.File(label="Reference FASTA or GenBank", file_count="single", type="filepath")
                    sample_files = gr.File(label="Sample files", file_count="multiple", type="filepath")
                with gr.Tab("Server paths"):
                    reference_path = gr.Textbox(
                        label="Reference path",
                        placeholder="/path/to/reference.fasta",
                    )
                    sample_paths = gr.Textbox(
                        label="Sample paths",
                        lines=6,
                        placeholder="/path/to/sample-directory\n/path/to/sample_R1.fastq.gz\n/path/to/sample_R2.fastq.gz",
                    )
        else:
            reference_file = gr.File(label="Reference FASTA or GenBank", file_count="single", type="filepath")
            sample_files = gr.File(label="Sample files", file_count="multiple", type="filepath")
            reference_path = gr.State(value="")
            sample_paths = gr.State(value="")

        with gr.Accordion("Advanced options", open=False):
            outdir_name = gr.Textbox(label="Output directory name", value=_default_outdir_name)
            temp_output_state = gr.State(value=bool(temp_output))
            prefix = gr.Textbox(label="Output prefix", value="snippy")
            cpus = gr.Number(label="Total CPUs", value=default_cpus, precision=0, minimum=1, maximum=max_cpus)
            cpus_per_sample = gr.Number(label="CPUs per sample", value=default_cpus_per_sample, precision=0, minimum=1)
            ram = gr.Number(label="RAM limit in GB", value=8, precision=0, minimum=1)
            max_depth = gr.Number(label="Maximum directory depth", value=4, precision=0, minimum=0)
            aggressive_ids = gr.Checkbox(label="Aggressive sample ID parsing", value=False)
            exclude_name_regex = gr.Textbox(label="Exclude filename regex", value=r"^(Undetermined|NTC|PTC)")
            core = gr.Slider(label="Core alignment threshold", minimum=0, maximum=1, value=0.95, step=0.01)
            inclusion_threshold = gr.Slider(label="Cluster inclusion threshold", minimum=0, maximum=1, value=0.1, step=0.01)
            stop_on_failure = gr.Checkbox(label="Stop on first sample failure", value=False)
            skip_check = gr.Checkbox(label="Skip dependency checks", value=False)
            force = gr.Checkbox(label="Force overwrite existing output directory", value=False)

        run_button = gr.Button("Run Snippy-NG", variant="primary")
        output_dir_state = gr.State(value=None)
        status = gr.Textbox(label="Run status", lines=5, visible=False)

        @gr.render(inputs=output_dir_state)
        def render_output_files(output_dir: str | None):
            if not output_dir:
                return
            tarball_button = gr.Button("Create tarball of output directory for download")
            tarball_file = gr.File(label="Download output tarball", visible=False, interactive=False)
            tarball_button.click(lambda: gr.update(visible=True), outputs=tarball_file, api_visibility="private")
            tarball_button.click(create_output_tarball, inputs=[output_dir_state, tarball_file], outputs=tarball_file).then(
                lambda: gr.update(visible=False), outputs=tarball_button, api_visibility="private"
            )

            file_explorer = gr.FileExplorer(
                label="Choose a file to download",
                file_count="single",
                root_dir=output_dir,
                glob="**/*",
                height=400,
                interactive=True,
            )
            output_file = gr.File(label="Download file", visible=False, interactive=False)
            file_explorer.change(lambda: gr.update(visible=True), outputs=output_file, api_visibility="private")
            file_explorer.change(select_download_file, file_explorer, output_file)

        run_event = run_button.click(
            _disable_run_button,
            outputs=run_button,
            api_visibility="private"
        ).then(
            lambda: gr.update(visible=True),
            outputs=status,
            api_visibility="private"
        )
        run_event.then(
            run_snippy_multi_gui,
            inputs=[
                reference_file,
                reference_path,
                sample_files,
                sample_paths,
                outdir_name,
                temp_output_state,
                prefix,
                cpus,
                cpus_per_sample,
                ram,
                max_depth,
                aggressive_ids,
                exclude_name_regex,
                core,
                inclusion_threshold,
                stop_on_failure,
                skip_check,
                force,
            ],
            outputs=[status, output_dir_state],
            api_name="snippy-ng",
        ).then(
            _enable_run_button,
            outputs=run_button,
            api_visibility="private"

        )

    return app
