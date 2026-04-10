"""Test gathering samples from multiple input directories."""
import gzip
from pathlib import Path
import pytest
from snippy_ng.utils.gather import gather_samples_config
from snippy_ng.utils.gather import guess_sample_id


def write_fastq(path: Path, seq: str = "ACGT", header: str = "@read") -> None:
    path.write_text(f"{header}\n{seq}\n+\nIIII\n")


def write_fastq_gz(path: Path, seq: str = "ACGT", header: str = "@read") -> None:
    with gzip.open(path, "wt") as fh:
        fh.write(f"{header}\n{seq}\n+\nIIII\n")


def assert_short_sample(samples, sample_id: str, *, left_suffix: str, right_suffix: str | None = None) -> None:
    assert sample_id in samples
    assert samples[sample_id]["type"] == "short"
    assert samples[sample_id]["left"].endswith(left_suffix)
    if right_suffix is None:
        assert samples[sample_id]["right"] is None
    else:
        assert samples[sample_id]["right"].endswith(right_suffix)


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
    samples = cfg["samples"]
    
    # Should have two different sample IDs
    assert cfg["reference"] is None
    assert len(samples) == 2
    assert "samples-JKD6159" in samples
    assert "data-JKD6159" in samples
    assert samples["samples-JKD6159"]["type"] == "asm"
    assert samples["data-JKD6159"]["type"] == "asm"


def test_gather_single_directory_no_prefix(tmp_path):
    """When scanning a single directory, don't prefix sample IDs."""
    samples_dir = tmp_path / "samples"
    samples_dir.mkdir()
    (samples_dir / "JKD6159.fasta").write_text(">contig1\nACGT\n")
    
    # Scan single directory - should not prefix
    cfg = gather_samples_config([samples_dir])
    samples = cfg["samples"]
    
    assert cfg["reference"] is None
    assert len(samples) == 1
    assert "JKD6159" in samples  # No prefix!
    assert samples["JKD6159"]["type"] == "asm"


def test_gather_multiple_directories_with_subdirs(tmp_path):
    """Files in duplicate subdirs should be disambiguated by parent path."""
    # Create samples/sample1/R1.fastq
    samples_dir = tmp_path / "samples"
    sample1_dir = samples_dir / "sample1"
    sample1_dir.mkdir(parents=True)
    write_fastq(sample1_dir / "R1.fastq")
    write_fastq(sample1_dir / "R2.fastq")
    
    # Create data/sample1/R1.fastq (same subdir name!)
    data_dir = tmp_path / "data"
    sample1b_dir = data_dir / "sample1"
    sample1b_dir.mkdir(parents=True)
    write_fastq(sample1b_dir / "R1.fastq")
    write_fastq(sample1b_dir / "R2.fastq")
    
    # Scan both - duplicate sample IDs should be disambiguated
    cfg = gather_samples_config([samples_dir, data_dir])
    samples = cfg["samples"]

    assert len(samples) == 2
    assert "samples-sample1" in samples
    assert "data-sample1" in samples
    assert samples["samples-sample1"]["type"] == "short"
    assert samples["data-sample1"]["type"] == "short"


def test_gather_same_filename_with_compression_keeps_extension(tmp_path):
    """Sample.fa and Sample.fa.gz should be treated as distinct sample IDs."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    (data_dir / "Sample.fa").write_text(">contig1\nACGT\n")
    with gzip.open(data_dir / "Sample.fa.gz", "wt") as fh:
        fh.write(">contig1\nACGT\n")

    cfg = gather_samples_config([data_dir])
    samples = cfg["samples"]

    assert "Sample.fa" in samples
    assert "Sample.fa.gz" in samples
    assert samples["Sample.fa"]["type"] == "asm"
    assert samples["Sample.fa.gz"]["type"] == "asm"


def test_gather_r1_r2_in_sample_dirs_uses_directory_name(tmp_path):
    """sample1/R1,R2 + sample2/R1,R2 should resolve to sample1 and sample2."""
    for sample in ("sample1", "sample2"):
        sample_dir = tmp_path / sample
        sample_dir.mkdir()
        write_fastq(sample_dir / "R1")
        write_fastq(sample_dir / "R2", seq="TGCA")

    cfg = gather_samples_config([tmp_path])
    samples = cfg["samples"]

    assert "sample1" in samples
    assert "sample2" in samples
    assert samples["sample1"]["type"] == "short"
    assert samples["sample2"]["type"] == "short"


def test_gather_mutant_r1_r2_in_sample_dirs_prefixes_parent(tmp_path):
    """sample1/mutant_R1,R2 + sample2/mutant_R1,R2 should be parent-disambiguated."""
    for sample in ("sample1", "sample2"):
        sample_dir = tmp_path / sample
        sample_dir.mkdir()
        write_fastq(sample_dir / "mutant_R1")
        write_fastq(sample_dir / "mutant_R2", seq="TGCA")

    cfg = gather_samples_config([tmp_path])
    samples = cfg["samples"]

    assert "sample1-mutant" in samples
    assert "sample2-mutant" in samples
    assert samples["sample1-mutant"]["type"] == "short"
    assert samples["sample2-mutant"]["type"] == "short"


def test_gather_same_filename_in_outbreak_dirs_prefixes_outbreak(tmp_path):
    """outbreak1/sample1 + outbreak2/sample1 should become outbreak1-sample1/outbreak2-sample1."""
    for outbreak in ("outbreak1", "outbreak2"):
        outbreak_dir = tmp_path / outbreak
        outbreak_dir.mkdir()
        (outbreak_dir / "sample1.fa").write_text(">contig1\nACGT\n")

    cfg = gather_samples_config([tmp_path])
    samples = cfg["samples"]

    assert "outbreak1-sample1" in samples
    assert "outbreak2-sample1" in samples
    assert samples["outbreak1-sample1"]["type"] == "asm"
    assert samples["outbreak2-sample1"]["type"] == "asm"


def test_gather_respects_symlink_paths(tmp_path):
    """Real and symlinked input roots should stay distinct in sample naming."""
    real_dir = tmp_path / "outbreak_real"
    real_dir.mkdir()
    (real_dir / "sample1.fa").write_text(">contig1\nACGT\n")

    link_dir = tmp_path / "outbreak_link"
    link_dir.symlink_to(real_dir, target_is_directory=True)

    cfg = gather_samples_config([real_dir, link_dir])
    samples = cfg["samples"]

    assert "outbreak_real-sample1" in samples
    assert "outbreak_link-sample1" in samples
    assert samples["outbreak_real-sample1"]["type"] == "asm"
    assert samples["outbreak_link-sample1"]["type"] == "asm"


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
    samples = cfg["samples"]

    assert "JKD6159.fastq.gz" in samples
    assert "JKD6159.fasta" in samples
    assert samples["JKD6159.fastq.gz"]["type"] == "long"
    assert samples["JKD6159.fasta"]["type"] == "asm"


def test_gather_reference_special_case_excluded_from_samples(tmp_path):
    """Reference should be returned separately and excluded from discovered samples."""
    ref = tmp_path / "reference.fa"
    ref.write_text(">ref\nACGT\n")

    # sample with a potentially clashing ID should still be kept as a sample
    sample_dir = tmp_path / "reads"
    sample_dir.mkdir()
    write_fastq(sample_dir / "reference.fastq")

    cfg = gather_samples_config([tmp_path], reference=ref)
    samples = cfg["samples"]

    assert cfg["reference"] == str(ref.absolute())
    assert "reads-reference" in samples
    assert "reference.fa" not in samples


def test_gather_reference_id_conflict_is_disambiguated(tmp_path):
    """If a sample ID clashes with reference ID, sample must be renamed, not dropped."""
    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref\nACGT\n")

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "ref.fasta").write_text(">contig\nACGT\n")

    cfg = gather_samples_config([data_dir], reference=ref)
    samples = cfg["samples"]

    assert cfg["reference"] == str(ref.absolute())
    assert "ref" not in samples
    assert "data-ref" in samples
    assert samples["data-ref"]["type"] == "asm"


def test_gather_trimmed_illumina_pair_groups_into_single_sample(tmp_path):
    """mutant_1.trim/mutant_2.trim FASTQs should be grouped as one short-read sample."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    for name in ("mutant_1.trim.fastq.gz", "mutant_2.trim.fastq.gz"):
        write_fastq_gz(data_dir / name)

    cfg = gather_samples_config([data_dir])
    samples = cfg["samples"]

    assert_short_sample(
        samples,
        "mutant",
        left_suffix="mutant_1.trim.fastq.gz",
        right_suffix="mutant_2.trim.fastq.gz",
    )
    assert samples["mutant"]["left"] == str(data_dir / "mutant_1.trim.fastq.gz")
    assert samples["mutant"]["right"] == str(data_dir / "mutant_2.trim.fastq.gz")
    assert "mutant_1.trim" not in samples
    assert "mutant_2.trim" not in samples


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("sample_1.fastq.gz", "sample"),
        ("sample_2.fastq.gz", "sample"),
        ("sample_R1.fastq.gz", "sample"),
        ("sample_R2.fastq.gz", "sample"),
        ("sample_R1_001.fastq.gz", "sample"),
        ("sample_1.trim.fastq.gz", "sample"),
        ("sample_2.trimmed.fastq.gz", "sample"),
        ("sample_R1.clean.fastq.gz", "sample"),
    ],
)
def test_guess_sample_id_strips_common_read_suffix_variants(filename, expected):
    assert guess_sample_id(filename) == expected


def test_gather_single_short_read_keeps_right_as_none(tmp_path):
    """A single short-read file should not serialize missing R2 as the string 'None'."""
    reads = tmp_path / "single.fastq"
    write_fastq(reads)

    cfg = gather_samples_config([tmp_path])
    samples = cfg["samples"]

    assert_short_sample(samples, "single", left_suffix="single.fastq")
    assert samples["single"]["left"] == str(reads)


def test_gather_srr_pair_with_numeric_sample_id_groups_r1_r2_correctly(tmp_path):
    """SRR-style names ending in _1/_2 should not be misclassified as duplicate R1s."""
    r1 = tmp_path / "SRR6171131_1.fastq"
    r2 = tmp_path / "SRR6171131_2.fastq"
    write_fastq(r1)
    write_fastq(r2, seq="TGCA")

    cfg = gather_samples_config([tmp_path])
    samples = cfg["samples"]

    assert_short_sample(samples, "SRR6171131", left_suffix="SRR6171131_1.fastq", right_suffix="SRR6171131_2.fastq")
    assert samples["SRR6171131"]["left"] == str(r1)
    assert samples["SRR6171131"]["right"] == str(r2)


def test_gather_kitchen_sink_mixed_layouts_and_naming_styles(tmp_path):
    """Gather should handle a realistic mix of roots, file kinds, pair styles, and collisions."""
    ref = tmp_path / "reference.fa"
    ref.write_text(">ref\nACGT\n")

    outbreak_a = tmp_path / "outbreak_a"
    outbreak_b = tmp_path / "outbreak_b"
    misc = tmp_path / "misc"
    outbreak_a.mkdir()
    outbreak_b.mkdir()
    misc.mkdir()

    # Paired short reads with explicit R1/R2 suffixes.
    sample_alpha = outbreak_a / "sample_alpha"
    sample_alpha.mkdir()
    write_fastq(sample_alpha / "alpha_R1.fastq")
    write_fastq(sample_alpha / "alpha_R2.fastq", seq="TGCA")

    # Paired short reads using trimmed _1/_2 naming.
    trimmed = outbreak_a / "trimmed"
    trimmed.mkdir()
    write_fastq_gz(trimmed / "mutant_1.trim.fastq.gz")
    write_fastq_gz(trimmed / "mutant_2.trimmed.fastq.gz", seq="TGCA")

    # SRR-style numeric sample id pair.
    write_fastq_gz(outbreak_a / "SRR6171131_1.fastq.gz")
    write_fastq_gz(outbreak_a / "SRR6171131_2.fastq.gz", seq="TGCA")

    # Directory-based sample naming with plain R1/R2 filenames.
    sample_1 = outbreak_b / "sample_1"
    sample_2 = outbreak_b / "sample_2"
    sample_1.mkdir()
    sample_2.mkdir()
    write_fastq_gz(sample_1 / "R1.fastq.gz")
    write_fastq_gz(sample_1 / "R2.fastq.gz", seq="TGCA")
    write_fastq_gz(sample_2 / "R1.fastq.gz")
    write_fastq_gz(sample_2 / "R2.fastq.gz", seq="TGCA")

    # Single short read should stay unpaired.
    write_fastq(misc / "single.clean.fastq")

    # ONT-like long read.
    write_fastq_gz(misc / "nanopore.fastq.gz", header="@550e8400-e29b-41d4-a716-446655440000")

    # Assemblies with same basename across roots should be prefixed.
    (outbreak_a / "sample1.fasta").write_text(">contig1\nACGT\n")
    (outbreak_b / "sample1.fasta").write_text(">contig2\nTGCA\n")

    # Sample that would clash with reference id should be disambiguated.
    write_fastq(outbreak_b / "reference.fastq")

    cfg = gather_samples_config([outbreak_a, outbreak_b, misc], reference=ref)
    samples = cfg["samples"]

    assert cfg["reference"] == str(ref.absolute())

    assert_short_sample(samples, "alpha", left_suffix="alpha_R1.fastq", right_suffix="alpha_R2.fastq")

    assert_short_sample(samples, "mutant", left_suffix="mutant_1.trim.fastq.gz", right_suffix="mutant_2.trimmed.fastq.gz")

    assert_short_sample(samples, "SRR6171131", left_suffix="SRR6171131_1.fastq.gz", right_suffix="SRR6171131_2.fastq.gz")

    assert_short_sample(samples, "sample_1", left_suffix="sample_1/R1.fastq.gz", right_suffix="sample_1/R2.fastq.gz")

    assert_short_sample(samples, "sample_2", left_suffix="sample_2/R1.fastq.gz", right_suffix="sample_2/R2.fastq.gz")

    assert_short_sample(samples, "single.clean", left_suffix="single.clean.fastq")

    assert "nanopore" in samples
    assert samples["nanopore"]["type"] == "long"
    assert samples["nanopore"]["reads"].endswith("nanopore.fastq.gz")

    assert "outbreak_a-sample1" in samples
    assert "outbreak_b-sample1" in samples
    assert samples["outbreak_a-sample1"]["type"] == "asm"
    assert samples["outbreak_b-sample1"]["type"] == "asm"

    assert "outbreak_b-reference" in samples
    assert samples["outbreak_b-reference"]["type"] == "short"
    assert "reference" not in samples
