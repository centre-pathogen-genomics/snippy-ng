import gzip


def guess_format(fname):
    try:
        # Try to open as text, if fails, try gzip
        try:
            with open(fname, 'rt') as fh:
                line = fh.readline()
        except UnicodeDecodeError:
            with gzip.open(fname, 'rt') as fh:
                line = fh.readline()
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
    except IOError:
        return None
