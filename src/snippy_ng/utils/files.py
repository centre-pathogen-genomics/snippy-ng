from pathlib import Path


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
