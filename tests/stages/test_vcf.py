from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.vcf import MIXED_SITE_GT_FILTER
from snippy_ng.stages.vcf import VariantContextFilter
from snippy_ng.stages.vcf import CollapseDiploidGenotypes
from snippy_ng.stages.vcf import VcfFilterLong
from snippy_ng.stages.vcf import VcfFilterShort
from snippy_ng.stages.vcf import VcfToTab


def _records(vcf: Path) -> list[list[str]]:
    return [
        line.strip().split("\t")
        for line in vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def test_variant_context_filter_marks_local_snp_clusters_lowqual(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t2\t.\tC\tT\t60\t.\t.",
                "chr1\t6\t.\tC\tT\t60\t.\t.",
                "chr1\t10\t.\tC\tT\t60\t.\t.",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "output.vcf"

    VariantContextFilter.filter_vcf(
        input_vcf,
        output_vcf,
        max_local_snps=1,
        local_snp_window=5,
    )

    records = _records(output_vcf)
    assert len(records) == 3
    assert all(record[6] == "LowQual" for record in records)
    assert all("CONTEXT_LOWQUAL_REASON=LOCAL_SNP_CLUSTER" in record[7] for record in records)


def test_variant_context_filter_marks_snps_near_indels_lowqual(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t10\t.\tA\tT\t60\t.\t.",
                "chr1\t15\t.\tAC\tA\t60\t.\t.",
                "chr1\t50\t.\tG\tC\t60\t.\t.",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "output.vcf"

    VariantContextFilter.filter_vcf(
        input_vcf,
        output_vcf,
        min_snp_distance_to_indel=10,
    )

    records = _records(output_vcf)
    assert [record[:5] for record in records] == [
        ["chr1", "10", ".", "A", "T"],
        ["chr1", "15", ".", "AC", "A"],
        ["chr1", "50", ".", "G", "C"],
    ]
    assert records[0][6] == "LowQual"
    assert "CONTEXT_LOWQUAL_REASON=NEAR_INDEL" in records[0][7]
    assert records[1][6] == "."
    assert records[2][6] == "."


def test_variant_context_filter_marks_snps_near_alignment_edges_lowqual(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t10\t.\tA\tT\t60\t.\tMUMMER_EDGE_DIST=1",
                "chr1\t50\t.\tG\tC\t60\t.\tMUMMER_EDGE_DIST=25",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "output.vcf"

    VariantContextFilter.filter_vcf(
        input_vcf,
        output_vcf,
        min_snp_distance_to_breakpoint=10,
    )

    records = _records(output_vcf)
    assert [record[:5] for record in records] == [
        ["chr1", "10", ".", "A", "T"],
        ["chr1", "50", ".", "G", "C"],
    ]
    assert records[0][6] == "LowQual"
    assert "CONTEXT_LOWQUAL_REASON=NEAR_ALIGNMENT_EDGE" in records[0][7]
    assert records[1][6] == "."


def test_variant_context_filter_does_not_duplicate_existing_headers(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##FILTER=<ID=LowQual,Description="Existing caller filter">',
                '##INFO=<ID=CONTEXT_LOWQUAL_REASON,Number=.,Type=String,Description="Existing context reasons">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t10\t.\tA\tT\t60\t.\t.",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "output.vcf"

    VariantContextFilter.filter_vcf(input_vcf, output_vcf)

    headers = [line for line in output_vcf.read_text(encoding="utf-8").splitlines() if line.startswith("##")]
    assert sum(line.startswith("##FILTER=<ID=LowQual,") for line in headers) == 1
    assert sum(line.startswith("##INFO=<ID=CONTEXT_LOWQUAL_REASON,") for line in headers) == 1


def test_variant_context_filter_does_nothing_when_all_filters_are_disabled(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t10\t.\tA\tT\t60\t.\t.",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    output_vcf = tmp_path / "output.vcf"
    stage = VariantContextFilter(
        vcf=input_vcf,
        max_local_snps=0,
        local_snp_window=0,
        min_snp_distance_to_indel=0,
        min_snp_distance_to_breakpoint=0,
    )

    assert stage.enabled is False
    stage.filter_vcf(input_vcf, output_vcf)
    assert output_vcf.read_text(encoding="utf-8") == input_vcf.read_text(encoding="utf-8")


def test_variant_context_filter_is_enabled_when_a_rule_is_configured(tmp_path):
    stage = VariantContextFilter(
        vcf=tmp_path / "input.vcf",
        min_snp_distance_to_indel=10,
    )

    assert stage.enabled is True


def test_vcf_to_tab_uses_bcftools_query_with_simple_default_columns():
    stage = VcfToTab(
        vcf=Path("sample.vcf"),
        prefix="snps",
    )

    commands = stage.create_commands(Context())

    assert len(commands) == 1
    assert commands[0].processes[0].command == [
        "bcftools",
        "query",
        "--allow-undef-tags",
        "-f",
        "%CHROM\\t%POS\\t%TYPE\\t%REF\\t%ALT\\t%INFO/BCSQ\\n",
        "sample.vcf",
    ]
    assert commands[0].processes[1].command[0] == "awk"
    assert 'print "CHROM","POS","TYPE","REF","ALT","Consequence","gene","transcript","biotype","strand","amino_acid_change","dna_change"' in commands[0].processes[1].command[1]
    assert 'split(selected, bcsq, "\\|")' in commands[0].processes[1].command[1]
    assert commands[0].output_file == Path("snps.tab")


def test_vcf_filter_short_marks_heterozygous_sites_as_mixedsite():
    stage = VcfFilterShort(
        vcf=Path("calls.vcf"),
        reference=Path("ref.fa"),
        prefix="snps",
    )

    commands = stage.create_commands(Context(ram=4))

    mixed_site_cmd = commands[0].processes[-1]
    assert mixed_site_cmd.command == [
        "bcftools",
        "filter",
        "-s",
        "MixedSite",
        "-m",
        "+",
        "-e",
        MIXED_SITE_GT_FILTER,
        "-",
    ]


def test_vcf_filter_short_splits_multiallelic_sites_before_collapse():
    stage = VcfFilterShort(
        vcf=Path("calls.vcf"),
        reference=Path("ref.fa"),
        prefix="snps",
    )

    commands = stage.create_commands(Context(ram=4))

    norm_cmd = commands[0].processes[1]
    assert norm_cmd.command == [
        "bcftools",
        "norm",
        "-f",
        "ref.fa",
        "--check-ref",
        "e",
        "--multiallelics",
        "-",
        "-Ou",
    ]


def test_vcf_filter_long_marks_heterozygous_sites_as_mixedsite():
    stage = VcfFilterLong(
        vcf=Path("calls.vcf"),
        reference=Path("ref.fa"),
        reference_index=Path("ref.fa.fai"),
        prefix="snps",
    )

    commands = stage.create_commands(Context(ram=4))

    mixed_site_cmd = commands[2].processes[-1]
    assert mixed_site_cmd.command == [
        "bcftools",
        "filter",
        "-s",
        "MixedSite",
        "-m",
        "+",
        "-e",
        MIXED_SITE_GT_FILTER,
        "-",
    ]


def test_vcf_filter_long_splits_multiallelic_sites_before_collapse():
    stage = VcfFilterLong(
        vcf=Path("calls.vcf"),
        reference=Path("ref.fa"),
        reference_index=Path("ref.fa.fai"),
        prefix="snps",
    )

    commands = stage.create_commands(Context(ram=4))

    norm_cmd = commands[2].processes[2]
    assert norm_cmd.command == [
        "bcftools",
        "norm",
        "-f",
        "ref.fa",
        "--check-ref",
        "e",
        "--multiallelics",
        "-",
        "-Ou",
    ]


def test_collapse_diploid_genotypes_rewrites_gt_field(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"
    input_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                "chr1\t10\t.\tA\tG\t60\tPASS\t.\tGT:DP\t1/1:12",
                "chr1\t11\t.\tA\tG\t60\tPASS\t.\tGT:DP\t0/1:9",
                "chr1\t12\t.\tA\tG\t60\tPASS\t.\tGT:DP\t0/0:10",
                "chr1\t13\t.\tA\tG\t60\tPASS\t.\tGT:DP\t1|1:7",
                "chr1\t14\t.\tA\tG\t60\tPASS\t.\tGT:DP\t1|0:5",
                "chr1\t15\t.\tA\tG\t60\tPASS\t.\tGT:DP\t0|0:3",
            ]
        ) + "\n",
        encoding="utf-8",
    )

    CollapseDiploidGenotypes.collapse_genotypes(input_vcf, output_vcf)

    records = _records(output_vcf)
    assert records[0][9] == "1:12"
    assert records[1][9] == ".:9"
    assert records[2][9] == "0:10"
    assert records[3][9] == "1:7"
    assert records[4][9] == ".:5"
    assert records[5][9] == "0:3"
