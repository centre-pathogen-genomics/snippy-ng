"""Tests for seq_utils module."""
import gzip
from snippy_ng.utils.seq import gather_samples_config, guess_reference_format


def test_guess_format_fasta(tmp_path):
    """Test detection of FASTA format."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">sequence1\nATGCATGC\n")
    assert guess_reference_format(str(fasta_file)) == "fasta"


def test_guess_format_genbank(tmp_path):
    """Test detection of GenBank format."""
    genbank_file = tmp_path / "test.gbk"
    genbank_file.write_text("LOCUS       AB000123\n")
    assert guess_reference_format(str(genbank_file)) == "genbank"


def test_guess_format_embl(tmp_path):
    """Test detection of EMBL format."""
    embl_file = tmp_path / "test.embl"
    embl_file.write_text("ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.\n")
    assert guess_reference_format(str(embl_file)) == "embl"


def test_guess_format_unknown(tmp_path):
    """Test detection returns None for unknown formats."""
    unknown_file = tmp_path / "test.txt"
    unknown_file.write_text("This is just some random text\n")
    assert guess_reference_format(str(unknown_file)) is None


def test_guess_format_empty_file(tmp_path):
    """Test detection returns None for empty files."""
    empty_file = tmp_path / "empty.txt"
    empty_file.write_text("")
    assert guess_reference_format(str(empty_file)) is None


def test_guess_format_gzipped_fasta(tmp_path):
    """Test detection of gzipped FASTA format."""
    gz_file = tmp_path / "test.fasta.gz"
    with gzip.open(gz_file, 'wt') as f:
        f.write(">sequence1\nATGCATGC\n")
    assert guess_reference_format(str(gz_file)) == "fasta"


def test_guess_format_gzipped_genbank(tmp_path):
    """Test detection of gzipped GenBank format."""
    gz_file = tmp_path / "test.gbk.gz"
    with gzip.open(gz_file, 'wt') as f:
        f.write("LOCUS       AB000123\n")
    assert guess_reference_format(str(gz_file)) == "genbank"


def test_guess_format_gzipped_embl(tmp_path):
    """Test detection of gzipped EMBL format."""
    gz_file = tmp_path / "test.embl.gz"
    with gzip.open(gz_file, 'wt') as f:
        f.write("ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.\n")
    assert guess_reference_format(str(gz_file)) == "embl"


def test_guess_format_nonexistent_file(tmp_path):
    """Test detection raises error for non-existent files."""
    nonexistent = tmp_path / "does_not_exist.txt"
    try:
        guess_reference_format(str(nonexistent))
    except FileNotFoundError:
        pass
    else:
        assert False, "Expected FileNotFoundError"


def test_gather_samples_config_ill_paired(tmp_path):
    """Detect Illumina paired reads and build config."""
    r1 = tmp_path / "sampleA_R1.fastq"
    r2 = tmp_path / "sampleA_R2.fastq"
    r1.write_text("@read1\nACGT\n+\n!!!!\n")
    r2.write_text("@read2\nTGCA\n+\n!!!!\n")

    cfg = gather_samples_config([tmp_path])
    assert "sampleA" in cfg
    assert cfg["sampleA"]["type"] == "short"
    assert cfg["sampleA"]["left"].endswith("sampleA_R1.fastq")
    assert cfg["sampleA"]["right"].endswith("sampleA_R2.fastq")


def test_gather_samples_config_ont_single(tmp_path):
    """Detect ONT reads by UUID-like header."""
    ont = tmp_path / "sampleB.fastq"
    ont.write_text("@550e8400-e29b-41d4-a716-446655440000\nACGT\n+\n!!!!\n")

    cfg = gather_samples_config([tmp_path])
    assert "sampleB" in cfg
    assert cfg["sampleB"]["type"] == "long"
    assert cfg["sampleB"]["reads"].endswith("sampleB.fastq")


def test_gather_samples_config_asm_single(tmp_path):
    """Detect assembly FASTA and build config."""
    asm = tmp_path / "sampleC.fasta"
    asm.write_text(">contig1\nACGT\n")

    cfg = gather_samples_config([tmp_path])
    assert "sampleC" in cfg
    assert cfg["sampleC"]["type"] == "asm"
    assert cfg["sampleC"]["assembly"].endswith("sampleC.fasta")


def test_gather_samples_config_excludes_by_name(tmp_path):
    """Exclude files by name regex."""
    r1 = tmp_path / "Undetermined_S0_R1.fastq"
    r2 = tmp_path / "Undetermined_S0_R2.fastq"
    r1.write_text("@read1\nACGT\n+\n!!!!\n")
    r2.write_text("@read2\nTGCA\n+\n!!!!\n")

    cfg = gather_samples_config([tmp_path])
    assert cfg == {}


def test_gather_samples_config_aggressive_ids(tmp_path):
    """Aggressive ID parsing should drop lane/sample suffixes."""
    r1 = tmp_path / "samp_L001_S1_R1.fastq"
    r2 = tmp_path / "samp_L001_S1_R2.fastq"
    r1.write_text("@read1\nACGT\n+\n!!!!\n")
    r2.write_text("@read2\nTGCA\n+\n!!!!\n")

    cfg = gather_samples_config([tmp_path], aggressive_ids=True)
    assert "samp" in cfg
    assert cfg["samp"]["type"] == "short"
