import gzip
from pathlib import Path


def _guess_reference_format_from_suffix(fname):
    suffixes = [suffix.lower() for suffix in Path(fname).suffixes]
    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]
    suffix = suffixes[-1] if suffixes else ""
    if suffix in {".fa", ".fna", ".fasta"}:
        return "fasta"
    if suffix in {".gb", ".gbk", ".gbff", ".genbank"}:
        return "genbank"
    if suffix == ".embl":
        return "embl"
    return None

def guess_reference_format(fname):
    # Try to open as text, if fails, try gzip
    try:
        fh = open(fname, 'rt')
        line = fh.readline()
        fh.close()
    except FileNotFoundError:
        return _guess_reference_format_from_suffix(fname)
    except UnicodeDecodeError:
        fh = gzip.open(fname, 'rt')
        line = fh.readline()
        fh.close()
    if not line:
        return None
    if line.startswith("LOCUS"):
        return 'genbank'
    elif line.startswith("ID "):
        return 'embl'
    elif line.startswith(">"):
        return 'fasta'
    else:
        return _guess_reference_format_from_suffix(fname)
