from pathlib import Path
import os

import pytest

from snippy_ng.context import Context
from snippy_ng.stages.calling import (
    Clair3Caller,
    Clair3ModelSelectorError,
    DoradoPolishCaller,
    FreebayesCaller,
    FreebayesCallerLong,
    LongbowClair3ModelSelector,
    MIN_SHORT_CHUNK_SIZE,
    PAFCaller,
    ShowSnpsCaller,
    get_short_chunk_size,
)


def test_get_short_chunk_size_scales_with_cpus_and_reference_index(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000000\t0\t0\t0\n")

    num_chunks, chunk_size = get_short_chunk_size(reference, reference_index, cpus=4)

    assert num_chunks == 7
    assert chunk_size == 142857


def test_get_short_chunk_size_honours_minimum(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000\t0\t0\t0\n")

    _, chunk_size = get_short_chunk_size(reference, reference_index, cpus=8)

    assert chunk_size == MIN_SHORT_CHUNK_SIZE


def test_freebayes_caller_uses_adaptive_chunk_size_for_region_generation(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000000\t0\t0\t0\n")
    bam = tmp_path / "reads.bam"
    bam.write_text("")
    bam_index = tmp_path / "reads.bam.bai"
    bam_index.write_text("")

    stage = FreebayesCaller(
        reference=reference,
        reference_index=reference_index,
        bam=bam,
        bam_index=bam_index,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=4))

    generate_regions = commands[0].processes[0]
    assert generate_regions.command == [
        "fasta_generate_regions.py",
        str(reference_index),
        "142857",
    ]


def test_freebayes_caller_uses_configured_mapping_quality(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000000\t0\t0\t0\n")
    bam = tmp_path / "reads.bam"
    bam.write_text("")
    bam_index = tmp_path / "reads.bam.bai"
    bam_index.write_text("")

    stage = FreebayesCaller(
        reference=reference,
        reference_index=reference_index,
        bam=bam,
        bam_index=bam_index,
        min_mapping_quality=30,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=4))
    freebayes_command = commands[1].command

    assert freebayes_command[
        freebayes_command.index("--min-mapping-quality") + 1
    ] == "30"


def test_freebayes_long_caller_uses_configured_mapping_quality(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1000000\t0\t0\t0\n")
    bam = tmp_path / "reads.bam"
    bam.write_text("")
    bam_index = tmp_path / "reads.bam.bai"
    bam_index.write_text("")

    stage = FreebayesCallerLong(
        reference=reference,
        reference_index=reference_index,
        bam=bam,
        bam_index=bam_index,
        min_mapping_quality=30,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=4))
    freebayes_command = commands[1].command

    assert freebayes_command[
        freebayes_command.index("--min-mapping-quality") + 1
    ] == "30"


def test_dorado_polish_caller_uses_bacteria_vcf_mode(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1\t0\t1\t2\n")
    bam = tmp_path / "reads.bam"
    bam.write_text("")
    bam_index = tmp_path / "reads.bam.bai"
    bam_index.write_text("")

    stage = DoradoPolishCaller(
        reference=reference,
        reference_index=reference_index,
        bam=bam,
        bam_index=bam_index,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=4))

    assert commands[0].command == [
        "dorado",
        "polish",
        str(bam),
        str(reference),
        "--vcf",
        "--threads",
        "4",
        "--bacteria",
    ]
    assert commands[0].output_file == Path("snippy.raw.vcf")


def test_paf_caller_disables_paftools_mapping_quality_prefilter(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1\t0\t1\t2\n")
    ref_dict = tmp_path / "ref.dict"
    ref_dict.write_text("chr1\t1\n")
    paf = tmp_path / "alignments.paf"
    paf.write_text("")

    stage = PAFCaller(
        paf=paf,
        ref_dict=ref_dict,
        reference=reference,
        reference_index=reference_index,
        min_mapping_quality=60,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=1))
    paftools_command = commands[2].command

    assert paftools_command[:2] == ["paftools.js", "call"]
    assert paftools_command[
        paftools_command.index("-q") + 1
    ] == "0"


def test_paf_caller_sorts_aligned_bed_by_reference_dict_before_complement(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1\t0\t1\t2\n")
    ref_dict = tmp_path / "ref.dict"
    ref_dict.write_text("chr1\t1\n")
    paf = tmp_path / "alignments.paf"
    paf.write_text("")

    stage = PAFCaller(
        paf=paf,
        ref_dict=ref_dict,
        reference=reference,
        reference_index=reference_index,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=1))
    missing_bed_pipeline = commands[1]

    assert [process.command for process in missing_bed_pipeline.processes] == [
        ["bedtools", "sort", "-g", str(ref_dict), "-i", "snippy.aln.bed"],
        ["bedtools", "complement", "-g", str(ref_dict), "-i", "-"],
    ]
    assert missing_bed_pipeline.output_file == Path("snippy.missing.bed")


def test_show_snps_caller_sorts_aligned_bed_by_reference_dict_before_complement(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1\t0\t1\t2\n")
    ref_dict = tmp_path / "ref.dict"
    ref_dict.write_text("chr1\t1\n")
    delta = tmp_path / "alignments.delta"
    delta.write_text("")
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">chr1\nA\n")

    stage = ShowSnpsCaller(
        delta=delta,
        ref_dict=ref_dict,
        assembly=assembly,
        reference=reference,
        reference_index=reference_index,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=1))
    missing_bed_pipeline = commands[3]

    assert [process.command for process in missing_bed_pipeline.processes] == [
        ["bedtools", "sort", "-g", str(ref_dict), "-i", "snippy.aln.bed"],
        ["bedtools", "complement", "-g", str(ref_dict), "-i", "-"],
    ]
    assert missing_bed_pipeline.output_file == Path("snippy.missing.bed")


def test_paf_vcf_annotation_marks_repetitive_mapping_lowqual(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsnippy",
                "chr1\t51\t.\tA\tG\t60\t.\t.\tGT\t1/1",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    paf = tmp_path / "alignments.paf"
    paf.write_text(
        "qry1\t1000\t0\t1000\t+\tchr1\t1000\t0\t100\t90\t100\t60"
        "\ttp:A:P\ts1:i:100\ts2:i:90\tcm:i:20\tAS:i:180\tdv:f:0.01\tde:f:0.01\trl:i:100\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "output.vcf"

    PAFCaller.annotate_paf_vcf(
        input_vcf,
        paf,
        output_vcf,
        min_mapping_quality=30,
        max_secondary_to_primary_ratio=0.8,
        max_gap_compressed_divergence=0.05,
        max_sequence_divergence=0.05,
        max_repeat_query_fraction=0.5,
        min_chain_minimizers_per_kb=0,
    )

    record = [
        line.strip().split("\t")
        for line in output_vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ][0]

    assert record[6] == "LowQual"
    assert "ASM_MAPQ=60" in record[7]
    assert "ASM_S2_RATIO=0.9" in record[7]
    assert "ASM_EDGE_DIST=49" in record[7]
    assert "ASM_LOWQUAL_REASON=REPETITIVE_CHAIN" in record[7]


def test_longbow_resolver_picks_r10_sup_model_from_env_root(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected_model = models_root / "r1041_e82_400bps_sup_v520"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell_version":"R10","basecaller":"Dorado","major_version":"Dorado0","basecalling_mode":"SUP","dorado_model_version":"V5.0.0"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_uses_launch_dir_for_relative_env_root(tmp_path, monkeypatch):
    launch_dir = tmp_path / "launch"
    launch_dir.mkdir()
    models_root = launch_dir / "clair3_models"
    expected_model = models_root / "r1041_e82_400bps_sup_v520"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"Sample":"reads.fastq.gz","Flowcell":"R10","Software":"guppy","Version":"5or6","Mode":"SUP"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("PWD", str(launch_dir))
    monkeypatch.setenv("CLAIR3_MODELS", "./clair3_models")

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_picks_r9_guppy5_legacy_model(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected_model = models_root / "ont_guppy5"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell":"R9","basecaller":"Guppy","major_version":"Guppy5/6","mode":"HAC"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_parses_title_case_longbow_json(tmp_path, monkeypatch):
    models_root = tmp_path / "models"
    expected_model = models_root / "r1041_e82_400bps_sup_v520"
    expected_model.mkdir(parents=True)

    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"Sample":"reads.fastq.gz","Flowcell":"R10","Software":"guppy","Version":"5or6","Mode":"SUP","Confidence level":"very high"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.setenv("CLAIR3_MODELS", str(models_root))

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == str(expected_model.resolve())


def test_longbow_resolver_errors_when_no_model_root_exists(tmp_path, monkeypatch):
    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell":"R10","basecaller":"Dorado","major_version":"Dorado0","mode":"SUP"}'
    )

    monkeypatch.delenv("CLAIR3_MODELS", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_ROOT", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_PATH", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_CONTAINER_ROOT", raising=False)
    monkeypatch.delenv("CLAIR3_MODELS_CONTAINER_ROOT", raising=False)
    monkeypatch.delenv("CONDA_PREFIX", raising=False)

    with pytest.raises(Clair3ModelSelectorError, match="Could not find any Clair3 model roots"):
        LongbowClair3ModelSelector.resolve_clair3_model(
            prediction_json,
            tmp_path / "resolved_model",
        )


def test_longbow_resolver_allows_container_only_model_root(tmp_path, monkeypatch):
    prediction_json = tmp_path / "longbow.json"
    prediction_json.write_text(
        '{"flowcell_version":"R10","basecaller":"Dorado","major_version":"Dorado0","basecalling_mode":"SUP","dorado_model_version":"V5.0.0"}'
    )
    output_path = tmp_path / "resolved_model"

    monkeypatch.delenv("CLAIR3_MODELS", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_ROOT", raising=False)
    monkeypatch.delenv("CLAIR3_MODEL_PATH", raising=False)
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    monkeypatch.setenv("CLAIR3_MODEL_CONTAINER_ROOT", "/opt/models")

    LongbowClair3ModelSelector.resolve_clair3_model(prediction_json, output_path)

    assert output_path.read_text(encoding="utf-8").strip() == "/opt/models/r1041_e82_400bps_sup_v500"


def test_clair3_caller_reads_model_path_from_manifest_file(tmp_path):
    manifest = tmp_path / "snippy.clair3_model.txt"
    manifest.write_text("/opt/models/r1041_e82_400bps_sup_v520\n", encoding="utf-8")
    bam = tmp_path / "reads.bam"
    bam.write_text("")
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n")
    reference_index = tmp_path / "ref.fa.fai"
    reference_index.write_text("chr1\t1\t0\t1\t2\n")

    stage = Clair3Caller(
        bam=bam,
        reference=reference,
        reference_index=reference_index,
        clair3_model=manifest,
        prefix="snippy",
    )

    commands = stage.create_commands(Context(cpus=1))

    assert commands[0].command[1] == "--model_path=/opt/models/r1041_e82_400bps_sup_v520"
    assert "--chunk_size=10000" in commands[0].command


def test_delta_to_vcf_converts_snps_and_indels(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(
        ">chr_snp\nACGTACGT\n>chr_del\nACGTACGT\n>chr_ins\nACGTACGT\n",
        encoding="utf-8",
    )
    assembly = tmp_path / "qry.fa"
    assembly.write_text(
        ">qry_snp\nATGTACGT\n>qry_del\nACGCGT\n>qry_ins\nACGTACGGT\n",
        encoding="utf-8",
    )
    delta = tmp_path / "calls.delta"
    delta.write_text(
        "\n".join(
            [
                f"{reference} {assembly}",
                "NUCMER",
                ">chr_snp qry_snp 8 8",
                "1 8 1 8 0 0 0",
                "0",
                ">chr_del qry_del 8 6",
                "1 8 1 6 0 0 0",
                "4",
                "1",
                "0",
                ">chr_ins qry_ins 8 9",
                "1 8 1 9 0 0 0",
                "-7",
                "0",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    tmp_vcf = tmp_path / "raw.vcf"
    final_vcf = tmp_path / "final.vcf"

    ShowSnpsCaller.delta_to_vcf(delta, reference, assembly, tmp_vcf)
    ShowSnpsCaller.postprocess_delta_vcf(tmp_vcf, final_vcf, "sample1")

    records = [
        line.strip().split("\t")
        for line in final_vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]

    assert [record[:5] for record in records] == [
        ["chr_del", "3", ".", "GTA", "G"],
        ["chr_ins", "6", ".", "C", "CG"],
        ["chr_snp", "2", ".", "C", "T"],
    ]
    assert all("MUMMER_BLEN=" in record[7] for record in records)
    assert all(record[8] == "GT:DP:AO" for record in records)
    assert all(record[9] == "1/1:1:1" for record in records)


def test_delta_to_vcf_filters_low_identity_blocks(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr_keep\nACGTACGT\n>chr_drop\nACGTACGT\n", encoding="utf-8")
    assembly = tmp_path / "qry.fa"
    assembly.write_text(">qry_keep\nATGTACGT\n>qry_drop\nATGTACGT\n", encoding="utf-8")
    delta = tmp_path / "calls.delta"
    delta.write_text(
        "\n".join(
            [
                f"{reference} {assembly}",
                "NUCMER",
                ">chr_keep qry_keep 8 8",
                "1 8 1 8 0 0 0",
                "0",
                ">chr_drop qry_drop 8 8",
                "1 8 1 8 2 0 0",
                "0",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "raw.vcf"

    ShowSnpsCaller.delta_to_vcf(
        delta,
        reference,
        assembly,
        output_vcf,
        min_delta_identity=90.0,
    )

    records = [
        line.strip().split("\t")
        for line in output_vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]

    assert [record[0] for record in records] == ["chr_keep"]


def test_delta_to_vcf_filters_ambiguous_overlapping_blocks(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    assembly = tmp_path / "qry.fa"
    assembly.write_text(">qry1\nATGTACGT\n>qry2\nATGTACGT\n", encoding="utf-8")
    delta = tmp_path / "calls.delta"
    delta.write_text(
        "\n".join(
            [
                f"{reference} {assembly}",
                "NUCMER",
                ">chr1 qry1 8 8",
                "1 8 1 8 0 0 0",
                "0",
                ">chr1 qry2 8 8",
                "1 8 1 8 0 0 0",
                "0",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "raw.vcf"

    ShowSnpsCaller.delta_to_vcf(delta, reference, assembly, output_vcf)

    records = [
        line
        for line in output_vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]

    assert records == []
