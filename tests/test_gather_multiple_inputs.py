"""Test gathering samples from multiple input directories."""
import gzip
from snippy_ng.utils.gather import gather_samples_config


def test_gather_multiple_directories_with_same_basename(tmp_path):
    """When scanning multiple directories with files of the same name, prefix with dir name."""
    # Create samples/JKD6159.fasta
    samples_dir = tmp_path / "samples"
    samples_dir.mkdir()
    (samples_dir / "JKD6159.fasta").write_text(">contig1\nACGT\n")
    
    # Create data/JKD6159.fasta
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "JKD6159.fasta").write_text(">contig2\nTGCA\n")
    
    # Scan both directories - should create separate samples
    cfg = gather_samples_config([samples_dir, data_dir])
    
    # Should have two different sample IDs
    assert len(cfg) == 2
    assert "samples-JKD6159" in cfg
    assert "data-JKD6159" in cfg
    assert cfg["samples-JKD6159"]["type"] == "asm"
    assert cfg["data-JKD6159"]["type"] == "asm"


def test_gather_single_directory_no_prefix(tmp_path):
    """When scanning a single directory, don't prefix sample IDs."""
    samples_dir = tmp_path / "samples"
    samples_dir.mkdir()
    (samples_dir / "JKD6159.fasta").write_text(">contig1\nACGT\n")
    
    # Scan single directory - should not prefix
    cfg = gather_samples_config([samples_dir])
    
    assert len(cfg) == 1
    assert "JKD6159" in cfg  # No prefix!
    assert cfg["JKD6159"]["type"] == "asm"


def test_gather_multiple_directories_with_subdirs(tmp_path):
    """Files in duplicate subdirs should be disambiguated by parent path."""
    # Create samples/sample1/R1.fastq
    samples_dir = tmp_path / "samples"
    sample1_dir = samples_dir / "sample1"
    sample1_dir.mkdir(parents=True)
    (sample1_dir / "R1.fastq").write_text("@read\nACGT\n+\nIIII\n")
    (sample1_dir / "R2.fastq").write_text("@read\nACGT\n+\nIIII\n")
    
    # Create data/sample1/R1.fastq (same subdir name!)
    data_dir = tmp_path / "data"
    sample1b_dir = data_dir / "sample1"
    sample1b_dir.mkdir(parents=True)
    (sample1b_dir / "R1.fastq").write_text("@read\nACGT\n+\nIIII\n")
    (sample1b_dir / "R2.fastq").write_text("@read\nACGT\n+\nIIII\n")
    
    # Scan both - duplicate sample IDs should be disambiguated
    cfg = gather_samples_config([samples_dir, data_dir])

    assert len(cfg) == 2
    assert "samples-sample1" in cfg
    assert "data-sample1" in cfg
    assert cfg["samples-sample1"]["type"] == "short"
    assert cfg["data-sample1"]["type"] == "short"


def test_gather_same_filename_with_compression_keeps_extension(tmp_path):
    """Sample.fa and Sample.fa.gz should be treated as distinct sample IDs."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    (data_dir / "Sample.fa").write_text(">contig1\nACGT\n")
    with gzip.open(data_dir / "Sample.fa.gz", "wt") as fh:
        fh.write(">contig1\nACGT\n")

    cfg = gather_samples_config([data_dir])

    assert "Sample.fa" in cfg
    assert "Sample.fa.gz" in cfg
    assert cfg["Sample.fa"]["type"] == "asm"
    assert cfg["Sample.fa.gz"]["type"] == "asm"


def test_gather_r1_r2_in_sample_dirs_uses_directory_name(tmp_path):
    """sample1/R1,R2 + sample2/R1,R2 should resolve to sample1 and sample2."""
    for sample in ("sample1", "sample2"):
        sample_dir = tmp_path / sample
        sample_dir.mkdir()
        (sample_dir / "R1").write_text("@read\nACGT\n+\nIIII\n")
        (sample_dir / "R2").write_text("@read\nTGCA\n+\nIIII\n")

    cfg = gather_samples_config([tmp_path])

    assert "sample1" in cfg
    assert "sample2" in cfg
    assert cfg["sample1"]["type"] == "short"
    assert cfg["sample2"]["type"] == "short"


def test_gather_mutant_r1_r2_in_sample_dirs_prefixes_parent(tmp_path):
    """sample1/mutant_R1,R2 + sample2/mutant_R1,R2 should be parent-disambiguated."""
    for sample in ("sample1", "sample2"):
        sample_dir = tmp_path / sample
        sample_dir.mkdir()
        (sample_dir / "mutant_R1").write_text("@read\nACGT\n+\nIIII\n")
        (sample_dir / "mutant_R2").write_text("@read\nTGCA\n+\nIIII\n")

    cfg = gather_samples_config([tmp_path])

    assert "sample1-mutant" in cfg
    assert "sample2-mutant" in cfg
    assert cfg["sample1-mutant"]["type"] == "short"
    assert cfg["sample2-mutant"]["type"] == "short"


def test_gather_same_filename_in_outbreak_dirs_prefixes_outbreak(tmp_path):
    """outbreak1/sample1 + outbreak2/sample1 should become outbreak1-sample1/outbreak2-sample1."""
    for outbreak in ("outbreak1", "outbreak2"):
        outbreak_dir = tmp_path / outbreak
        outbreak_dir.mkdir()
        (outbreak_dir / "sample1.fa").write_text(">contig1\nACGT\n")

    cfg = gather_samples_config([tmp_path])

    assert "outbreak1-sample1" in cfg
    assert "outbreak2-sample1" in cfg
    assert cfg["outbreak1-sample1"]["type"] == "asm"
    assert cfg["outbreak2-sample1"]["type"] == "asm"


def test_gather_respects_symlink_paths(tmp_path):
    """Real and symlinked input roots should stay distinct in sample naming."""
    real_dir = tmp_path / "outbreak_real"
    real_dir.mkdir()
    (real_dir / "sample1.fa").write_text(">contig1\nACGT\n")

    link_dir = tmp_path / "outbreak_link"
    link_dir.symlink_to(real_dir, target_is_directory=True)

    cfg = gather_samples_config([real_dir, link_dir])

    assert "outbreak_real-sample1" in cfg
    assert "outbreak_link-sample1" in cfg
    assert cfg["outbreak_real-sample1"]["type"] == "asm"
    assert cfg["outbreak_link-sample1"]["type"] == "asm"


def test_gather_mixed_kinds_keep_extension_with_cross_parent_collision(tmp_path):
    """Keep filename extensions for mixed-kind clashes even when sid also exists elsewhere."""
    tests_data = tmp_path / "tests" / "data"
    tests_data.mkdir(parents=True)
    (tests_data / "JKD6159.fasta").write_text(">contig1\nACGT\n")
    with gzip.open(tests_data / "JKD6159.fastq.gz", "wt") as fh:
        fh.write("@550e8400-e29b-41d4-a716-446655440000\nACGT\n+\nIIII\n")

    long_ref = tmp_path / "long" / "reference"
    long_ref.mkdir(parents=True)
    (long_ref / "JKD6159.fa").write_text(">contig2\nTGCA\n")

    cfg = gather_samples_config([tmp_path])

    assert "JKD6159.fastq.gz" in cfg
    assert "JKD6159.fasta" in cfg
    assert cfg["JKD6159.fastq.gz"]["type"] == "long"
    assert cfg["JKD6159.fasta"]["type"] == "asm"
