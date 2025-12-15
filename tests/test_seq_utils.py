"""Tests for seq_utils module."""
import gzip
from snippy_ng.seq_utils import guess_format


def test_guess_format_fasta(tmp_path):
    """Test detection of FASTA format."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">sequence1\nATGCATGC\n")
    assert guess_format(str(fasta_file)) == "fasta"


def test_guess_format_genbank(tmp_path):
    """Test detection of GenBank format."""
    genbank_file = tmp_path / "test.gbk"
    genbank_file.write_text("LOCUS       AB000123\n")
    assert guess_format(str(genbank_file)) == "genbank"


def test_guess_format_embl(tmp_path):
    """Test detection of EMBL format."""
    embl_file = tmp_path / "test.embl"
    embl_file.write_text("ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.\n")
    assert guess_format(str(embl_file)) == "embl"


def test_guess_format_unknown(tmp_path):
    """Test detection returns None for unknown formats."""
    unknown_file = tmp_path / "test.txt"
    unknown_file.write_text("This is just some random text\n")
    assert guess_format(str(unknown_file)) is None


def test_guess_format_empty_file(tmp_path):
    """Test detection returns None for empty files."""
    empty_file = tmp_path / "empty.txt"
    empty_file.write_text("")
    assert guess_format(str(empty_file)) is None


def test_guess_format_gzipped_fasta(tmp_path):
    """Test detection of gzipped FASTA format."""
    gz_file = tmp_path / "test.fasta.gz"
    with gzip.open(gz_file, 'wt') as f:
        f.write(">sequence1\nATGCATGC\n")
    assert guess_format(str(gz_file)) == "fasta"


def test_guess_format_gzipped_genbank(tmp_path):
    """Test detection of gzipped GenBank format."""
    gz_file = tmp_path / "test.gbk.gz"
    with gzip.open(gz_file, 'wt') as f:
        f.write("LOCUS       AB000123\n")
    assert guess_format(str(gz_file)) == "genbank"


def test_guess_format_gzipped_embl(tmp_path):
    """Test detection of gzipped EMBL format."""
    gz_file = tmp_path / "test.embl.gz"
    with gzip.open(gz_file, 'wt') as f:
        f.write("ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.\n")
    assert guess_format(str(gz_file)) == "embl"


def test_guess_format_nonexistent_file(tmp_path):
    """Test detection returns None for non-existent files."""
    nonexistent = tmp_path / "does_not_exist.txt"
    assert guess_format(str(nonexistent)) is None
