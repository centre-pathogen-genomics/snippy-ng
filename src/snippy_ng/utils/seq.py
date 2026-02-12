import gzip

def guess_reference_format(fname):
    # Try to open as text, if fails, try gzip
    try:
        fh = open(fname, 'rt')
        line = fh.readline()
        fh.close()
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
        return None