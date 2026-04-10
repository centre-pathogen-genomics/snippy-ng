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


def inner_join_csv_tsv(left_file, right_file, output_file, on):
    """
    Inner join two CSV/TSV files on a column.

    Parameters
    ----------
    left_file : str
    right_file : str
    output_file : str
    on : str
        column name to join on
    """

    def delimiter(path):
        return "\t" if Path(path).suffix.lower() == ".tsv" else ","

    left_delim = delimiter(left_file)
    right_delim = delimiter(right_file)
    out_delim = delimiter(output_file)

    # Load right file into lookup table
    with open(right_file, newline="", encoding="utf-8") as f:
        right_reader = csv.DictReader(f, delimiter=right_delim)
        right_rows = list(right_reader)
        right_fields = right_reader.fieldnames

    lookup = {}
    for row in right_rows:
        lookup.setdefault(row[on], []).append(row)

    # Stream left file and write joined output
    with open(left_file, newline="", encoding="utf-8") as lf, \
         open(output_file, "w", newline="", encoding="utf-8") as out:

        left_reader = csv.DictReader(lf, delimiter=left_delim)
        left_fields = left_reader.fieldnames

        right_extra = [c for c in right_fields if c != on]

        writer = csv.DictWriter(
            out,
            fieldnames=left_fields + right_extra,
            delimiter=out_delim
        )
        writer.writeheader()

        for left_row in left_reader:
            key = left_row[on]
            if key in lookup:
                for r in lookup[key]:
                    merged = left_row.copy()
                    for c in right_extra:
                        merged[c] = r[c]
                    writer.writerow(merged)