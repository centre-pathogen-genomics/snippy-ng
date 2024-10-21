def guess_format(fname):
    try:
        with open(fname, 'r') as fh:
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
