"""Select a Clair3 model from ONT read headers.

Dorado and Guppy encode the basecalling configuration in the read header, so the
matching Clair3 model can usually be chosen without running Longbow.  The
selection mirrors ``artic.utils.choose_model`` but resolves model directories
from the local Clair3 model roots instead of a fixed manifest.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional

from snippy_ng.__about__ import __version__
from snippy_ng.exceptions import Clair3ModelDownloadError, Clair3ModelSelectionError
from snippy_ng.logging import logger

CLAIR3_MODEL_BASE_URL = "https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch"
CLAIR3_MODEL_FILES = ("pileup.pt", "full_alignment.pt")
_MODEL_NAME_RE = re.compile(r"^[A-Za-z0-9_.+-]+$")


def candidate_model_roots() -> list[Path]:
    """Return the local directories that may contain Clair3 model directories."""
    launch_dir = Path(os.environ.get("PWD", str(Path.cwd()))).expanduser()
    roots: list[Path] = []
    for env_var in ("CLAIR3_MODELS", "CLAIR3_MODEL_ROOT"):
        value = os.environ.get(env_var)
        if value:
            root = Path(value).expanduser()
            if not root.is_absolute():
                root = launch_dir / root
            roots.append(root)

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        roots.extend(
            [
                Path(conda_prefix) / "bin" / "models",
                Path(conda_prefix) / "models",
            ]
        )

    roots.extend(
        [
            Path("/opt/models"),
            Path("/opt/models/clair3_models"),
        ]
    )

    unique_roots: list[Path] = []
    seen: set[Path] = set()
    for root in roots:
        root = root.expanduser()
        if root in seen:
            continue
        seen.add(root)
        unique_roots.append(root)
    return unique_roots


def container_model_roots() -> list[Path]:
    """Return model roots that are only valid inside an execution container."""
    roots: list[Path] = []
    for env_var in ("CLAIR3_MODEL_CONTAINER_ROOT", "CLAIR3_MODELS_CONTAINER_ROOT"):
        value = os.environ.get(env_var)
        if value:
            roots.append(Path(value))
    return roots


def _iter_model_directories(root: Path):
    """Yield model directories directly under *root* or *root*/clair3_models."""
    for base in (root, root / "clair3_models"):
        if not base.is_dir():
            continue
        for child in sorted(base.iterdir()):
            if child.is_dir() and not child.name.startswith("."):
                yield child


def local_model_index(roots: Optional[list[Path]] = None) -> dict[str, Path]:
    """Return ``{model_name: path}`` for models reachable from the local roots."""
    index: dict[str, Path] = {}
    for root in roots if roots is not None else candidate_model_roots():
        if not root.is_dir():
            continue
        for model_dir in _iter_model_directories(root):
            index.setdefault(model_dir.name, model_dir.resolve())
    return index


def ensure_models_available(roots: Optional[list[Path]] = None) -> dict[str, Path]:
    """Return the local model index, raising if no models can be used at all."""
    models = local_model_index(roots)
    if models:
        return models

    if roots is None and container_model_roots():
        # Cannot enumerate container-only roots from here; the caller falls
        # back to a container-aware path.
        return models

    searched = roots if roots is not None else candidate_model_roots()
    if not searched:
        raise Clair3ModelSelectionError(
            "No Clair3 model roots found. Set CLAIR3_MODELS to a directory of "
            "Clair3 models or provide a model with --model."
        )

    raise Clair3ModelSelectionError(
        "No Clair3 models found in "
        + ", ".join(str(root) for root in searched)
        + ". Download a model into one of these directories, set CLAIR3_MODELS, "
        "or provide a model with --model."
    )


def get_basecall_model_id(description: str) -> Optional[str]:
    """Extract the basecall model string from a FASTQ read description.

    Supports both Dorado header formats:
    - key=value: ``basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0``
    - SAM tags:  ``RG:Z:<UUID>_dna_r10.4.1_e8.2_400bps_hac@v5.2.0_barcode01``
    """
    tokens = description.split()

    key_values = {key: value for key, value in (token.split("=", 1) for token in tokens if "=" in token)}
    if "basecall_model_version_id" in key_values:
        return key_values["basecall_model_version_id"]

    for token in tokens:
        if token.startswith("RG:Z:"):
            rg_value = token[5:]
            match = re.search(r"(?:dna|rna)\w*_[^@\s]+@v[\d.]+", rg_value)
            if match:
                return match.group(0)

    return None


def read_first_read_header(reads: Path) -> Optional[str]:
    """Return the first non-empty header line of a FASTQ/FASTA file."""
    opener = gzip.open if reads.suffix.lower() == ".gz" else open
    with opener(reads, "rt", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line:
                return line
    return None


def _parse_basecall_model_id(model_id: str) -> Optional[dict[str, Optional[str]]]:
    parts = model_id.split("_")
    if len(parts) not in (4, 5):
        return None

    pore_type = parts[1].replace(".", "").lower()
    if len(parts) == 5:
        kit_id = parts[2].replace(".", "").lower()
        speed = parts[3].lower()
        model = parts[4]
        search_string = f"{pore_type}_{kit_id}_{speed}"
    else:
        model = parts[3]
        search_string = pore_type

    preset: Optional[str] = None
    version: Optional[str] = None
    if "@" in model:
        preset, basecaller_version = model.split("@", 1)
        preset = preset.lower()
        version = re.sub(r"\D", "", basecaller_version) or None
    else:
        preset = model.lower()

    return {
        "search_string": search_string,
        "preset": preset,
        "version": version,
        "model_id": model_id,
    }


def _select_model_name(
    model_id: str,
    names: list[str],
    model_directories: list[Path],
) -> str:
    def selection_error(
        message: str, *, suggest_download: bool = True
    ) -> Clair3ModelSelectionError:
        locations = ", ".join(str(path) for path in model_directories) or "none"
        available = ", ".join(sorted(names)) or "none"
        download_hint = ""
        if suggest_download:
            try:
                model_name = basecall_model_id_to_model_name(model_id)
            except Clair3ModelDownloadError:
                pass
            else:
                download_hint = (
                    " Or download it with: "
                    f"snippy-ng utils model download {model_name}."
                )
        return Clair3ModelSelectionError(
            f"{message} Searched model directories: {locations}. "
            f"Available models: {available}. Provide a model with --model.{download_hint}"
        )

    parsed = _parse_basecall_model_id(model_id)
    if parsed is None:
        raise selection_error(
            f"Could not parse basecall model id '{model_id}'.",
            suggest_download=False,
        )

    search_string = str(parsed["search_string"])
    preset = parsed["preset"]
    version = parsed["version"]

    candidates = [name for name in names if search_string in name.lower()]
    if not candidates:
        raise selection_error(
            f"No Clair3 model found for basecall model id '{model_id}' in the models directory."
        )

    if preset:
        candidates = [name for name in candidates if preset in name.lower()]
        if not candidates:
            raise selection_error(
                f"No Clair3 model found for basecall model id '{model_id}' in the models directory. "
                f"No model matches preset '{preset}'."
            )

    if len(candidates) == 1:
        return candidates[0]

    if version:
        marker = f"v{version}"
        by_version = [name for name in candidates if marker in name.lower()]
        if not by_version:
            raise selection_error(
                f"No Clair3 model found for basecall model id '{model_id}' in the models directory. "
                f"Candidates {', '.join(candidates)} do not match basecaller version v{version}."
            )
        candidates = by_version

    if len(candidates) == 1:
        return candidates[0]

    raise selection_error(
        f"Multiple Clair3 models match basecall model id '{model_id}': "
        f"{', '.join(candidates)}.",
        suggest_download=False,
    )


def choose_model(reads: Path, roots: Optional[list[Path]] = None) -> Optional[Path]:
    """Choose a Clair3 model for *reads* from the local model roots.

    Returns ``None`` when the read header has no basecall model id, in which
    case the caller should fall back to Longbow.  Raises
    :class:`Clair3ModelSelectionError` when a model id is present but no
    matching model directory exists.
    """
    header = read_first_read_header(reads)
    model_id = get_basecall_model_id(header) if header else None
    if model_id is None:
        return None

    models = ensure_models_available(roots)
    if not models:
        # Only container-only roots are configured; they cannot be inspected
        # here, so let the container-aware Longbow path resolve the model.
        return None

    model_directories = sorted({path.parent for path in models.values()})
    model_name = _select_model_name(model_id, list(models), model_directories)
    return models[model_name]


def basecall_model_id_to_model_name(value: str) -> str:
    """Convert a basecall model id or model directory name to a model name."""
    value = value.strip()
    if not value:
        raise Clair3ModelDownloadError("Model name or basecall model id is empty.")

    if "@" not in value:
        if value in {".", ".."} or not _MODEL_NAME_RE.fullmatch(value):
            raise Clair3ModelDownloadError(
                f"'{value}' is not a valid Clair3 model directory name."
            )
        return value

    if len(value.split("_")) != 5:
        raise Clair3ModelDownloadError(
            f"Could not derive a Clair3 model name from '{value}'. "
            "Provide a model directory name, for example r1041_e82_400bps_hac_v500."
        )

    parsed = _parse_basecall_model_id(value)
    if parsed is None or not parsed["version"]:
        raise Clair3ModelDownloadError(
            f"Could not derive a Clair3 model name from '{value}'. "
            "Provide a model directory name, for example r1041_e82_400bps_hac_v500."
        )

    model_name = f"{parsed['search_string']}_{parsed['preset']}_v{parsed['version']}"
    if not _MODEL_NAME_RE.fullmatch(model_name):
        raise Clair3ModelDownloadError(
            f"Derived model name '{model_name}' is not a valid directory name."
        )
    return model_name


def download_clair3_model(
    model: str,
    output_directory: Path,
    *,
    force: bool = False,
) -> Path:
    """Download a Clair3 model into *output_directory* and return its directory."""
    model_name = basecall_model_id_to_model_name(model)
    output_directory = Path(output_directory).expanduser()
    target_dir = output_directory / model_name

    expected_files = [target_dir / name for name in CLAIR3_MODEL_FILES]
    if not force and all(path.is_file() and path.stat().st_size > 0 for path in expected_files):
        logger.info(f"Clair3 model '{model_name}' already exists at {target_dir}")
        return target_dir

    output_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=output_directory, prefix=f".{model_name}-"
    ) as tmpdir:
        tmp_path = Path(tmpdir)
        for filename in CLAIR3_MODEL_FILES:
            url = f"{CLAIR3_MODEL_BASE_URL}/{model_name}/{filename}"
            destination = tmp_path / filename
            request = urllib.request.Request(
                url, headers={"User-Agent": f"snippy-ng/{__version__}"}
            )
            try:
                with urllib.request.urlopen(request) as response, destination.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    raise Clair3ModelDownloadError(
                        f"Clair3 model '{model_name}' is not available at {CLAIR3_MODEL_BASE_URL}. "
                        "Check the model name (for example r1041_e82_400bps_hac_v500)."
                    ) from exc
                raise Clair3ModelDownloadError(
                    f"Could not download Clair3 model '{model_name}' from '{url}' "
                    f"(HTTP {exc.code})."
                ) from exc
            except urllib.error.URLError as exc:
                raise Clair3ModelDownloadError(
                    f"Could not reach '{url}': {exc.reason}"
                ) from exc

            if destination.stat().st_size == 0:
                raise Clair3ModelDownloadError(
                    f"Downloaded Clair3 model file '{filename}' is empty."
                )

        target_dir.mkdir(parents=True, exist_ok=True)
        for filename in CLAIR3_MODEL_FILES:
            shutil.move(str(tmp_path / filename), str(target_dir / filename))

    logger.info(f"Downloaded Clair3 model '{model_name}' to {target_dir}")
    return target_dir
