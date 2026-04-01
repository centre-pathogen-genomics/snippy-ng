from typing import Optional
from pathlib import Path
import json
import csv


def load_metadata_as_json_str(metadata_path: Path) -> Optional[str]:
    metadata = None
    if metadata_path.suffix.lower() == ".json":
    # if metadata is a JSON file, read it and pass the content as a string to the pipeline
        with open(metadata_path, "r") as f:
            metadata = json.dumps(json.load(f))
    elif metadata_path.suffix.lower() == ".csv" or metadata_path.suffix.lower() == ".tsv":
        # if metadata is a CSV/TSV file, read it and convert to a list of dicts, then pass as a JSON string to the pipeline
        with open(metadata_path, "r") as f:
            if metadata_path.suffix.lower() == ".csv":
                reader = csv.DictReader(f)
            else:
                reader = csv.DictReader(f, delimiter="\t")
            metadata = json.dumps(list(reader))
    else:
        raise ValueError("Metadata file must be in JSON, CSV, or TSV format")
    return metadata

def human_readable_size(path: Path) -> str:
    """Return file/directory size in human-readable form (base 1024)."""
    if not path.exists():
        return ""

    try:
        if path.is_dir():
            total_bytes = sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
        else:
            total_bytes = path.stat().st_size
    except OSError:
        return ""

    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    size = float(total_bytes)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            if unit == "B":
                return f"{int(size)} {unit}"
            return f"{size:.1f} {unit}"
        size /= 1024

    return ""
