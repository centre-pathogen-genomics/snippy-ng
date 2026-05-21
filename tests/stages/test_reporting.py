import json
import sys
from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.reporting import SampleReport, TreeReport


def test_epi_report_preserves_small_branch_lengths():
    report = TreeReport()
    context = {
        "NEWICK": "(sample_a:0.00000004,sample_b:0.1);",
        "REPORT_NAME": "test-report",
        "METADATA": json.dumps(
            [
                {"id": "sample_a"},
                {"id": "sample_b"},
            ]
        ),
        "LOGS": "",
    }

    report.validate_context(context)

    assert "sample_a:0.0000000400" in context["NEWICK"]


def test_tree_report_render_converts_metadata_path_during_validation(tmp_path):
    tree = tmp_path / "tree.nwk"
    metadata = tmp_path / "metadata.csv"
    logs = tmp_path / "LOG.txt"
    template = tmp_path / "template.html"

    tree.write_text("(sample_a:0.1,sample_b:0.2);\n")
    metadata.write_text("sample,group\nsample_a,A\nsample_b,B\n")
    logs.write_text("tree report log\n")
    template.write_text("{{NEWICK}}\n{{METADATA_JSON}}\n{{LOGS}}\n{{COLOR_BY_COLUMN}}")

    report = TreeReport(
        template_path=template,
        context={
            "NEWICK": tree,
            "REPORT_NAME": "test-report",
            "METADATA": metadata,
            "LOGS": logs,
            "COLOR_BY_COLUMN": "group",
        },
        prefix=str(tmp_path / "tree-report"),
    )

    report.render_template()

    rendered = (tmp_path / "tree-report.html").read_text()
    assert "sample_a" in rendered
    assert '"group": "A"' in rendered
    assert "tree report log" in rendered
    assert "group" in rendered


def test_sample_report_parses_vcf_records_and_scope(tmp_path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type">\n'
        '##INFO=<ID=BCSQ,Number=.,Type=String,Description="BCF consequence annotation">\n'
        '##INFO=<ID=SB_ALT_FWD_FRAC,Number=1,Type=Float,Description="Alt forward fraction">\n'
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Info depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Fraction">\n'
        '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate observations">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
        "chr1\t10\t.\tA\tG\t99\tPASS\tTYPE=SNP;BCSQ=missense|gene1|x;SB_ALT_FWD_FRAC=0.75\tGT:DP:AF\t1/1:20:0.9\n"
        "chr1\t20\t.\tAT\tA\t12\tLowQual\tDP=5\tGT:AO\t0/1:2\n"
    )

    pass_records, format_fields, info_fields = SampleReport.parse_vcf_records(vcf, "pass")
    all_records, all_format_fields, all_info_fields = SampleReport.parse_vcf_records(vcf, "all")

    assert len(pass_records) == 1
    assert len(all_records) == 2
    assert [field["id"] for field in format_fields] == ["GT", "DP", "AF"]
    assert [field["id"] for field in all_format_fields] == ["GT", "DP", "AF", "AO"]
    assert [field["id"] for field in info_fields] == ["TYPE", "BCSQ", "SB_ALT_FWD_FRAC"]
    assert [field["id"] for field in all_info_fields] == ["TYPE", "BCSQ", "SB_ALT_FWD_FRAC", "DP"]
    assert info_fields[2]["field"] == "INFO:SB_ALT_FWD_FRAC"
    assert info_fields[2]["type"] == "Float"
    assert format_fields[1]["type"] == "Integer"
    assert pass_records[0]["sample"] == "sample1"
    assert pass_records[0]["type"] == "SNP"
    assert pass_records[0]["GT"] == "1/1"
    assert pass_records[0]["DP"] == 20
    assert pass_records[0]["AF"] == 0.9
    assert "depth" not in pass_records[0]
    assert "allele_fraction" not in pass_records[0]
    assert "sb_alt_fwd_frac" not in pass_records[0]
    assert pass_records[0]["INFO:SB_ALT_FWD_FRAC"] == 0.75
    assert pass_records[0]["consequence"] == "missense"


def test_sample_report_writes_variant_json_and_bed_windows(tmp_path):
    vcf = tmp_path / "sample.vcf"
    variants_json = tmp_path / "variants.json"
    regions_bed = tmp_path / "regions.bed"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t3\t.\tAT\tA\t99\tPASS\t.\n"
    )

    SampleReport.write_variant_assets(
        vcf,
        variants_json,
        regions_bed,
        variant_scope="pass",
        window_size=10,
    )

    payload = json.loads(variants_json.read_text())
    assert payload["variants"][0]["pos"] == 3
    assert payload["sampleFields"] == []
    assert payload["infoFields"] == []
    assert regions_bed.read_text() == "chr1\t0\t15\n"


def test_sample_report_render_embeds_payloads(tmp_path):
    template = tmp_path / "template.html"
    output = tmp_path / "report.html"
    variants = tmp_path / "variants.json"
    vcf = tmp_path / "sample.vcf"
    reference = tmp_path / "ref.fa"
    reference_index = tmp_path / "ref.fa.fai"
    cram = tmp_path / "sample.cram"
    crai = tmp_path / "sample.cram.crai"

    template.write_text("{{REPORT_NAME}} {{VARIANTS_JSON_B64}} {{VCF_B64}} {{CRAM_B64}} {{CRAI_B64}} {{REFERENCE_FASTA_B64}} {{REFERENCE_INDEX_B64}} {{REFERENCE_NAME}} {{VCF_NAME}} {{CRAM_NAME}} {{CRAI_NAME}} {{SAMPLE_NAME}} {{DATETIME}} {{VERSION}}")
    variants.write_text("[]")
    vcf.write_text("##fileformat=VCFv4.2\n")
    reference.write_text(">chr1\nA\n")
    reference_index.write_text("chr1\t1\t6\t1\t2\n")
    cram.write_bytes(b"cram")
    crai.write_bytes(b"crai")

    SampleReport.render_sample_report(
        template,
        output,
        variants,
        cram,
        crai,
        vcf,
        reference,
        reference_index,
        "Report",
        "sample1",
    )

    html = output.read_text()
    assert "Report" in html
    assert "sample1" in html
    assert "Y3JhbQ==" in html
    assert "sample.vcf" in html


def test_sample_report_render_without_igv_assets(tmp_path):
    template = tmp_path / "template.html"
    output = tmp_path / "report.html"
    variants = tmp_path / "variants.json"
    vcf = tmp_path / "sample.vcf"

    template.write_text("{{HAS_IGV}} {{CRAM_B64}} {{REFERENCE_FASTA_B64}} {{VCF_NAME}}")
    variants.write_text("[]")
    vcf.write_text("##fileformat=VCFv4.2\n")

    SampleReport.render_sample_report(
        template,
        output,
        variants,
        None,
        None,
        vcf,
        None,
        None,
        "Report",
    )

    assert output.read_text() == "false   sample.vcf"


def test_sample_report_commands_window_and_index_cram(tmp_path):
    stage = SampleReport(
        vcf=tmp_path / "sample.vcf",
        alignment=tmp_path / "sample.bam",
        reference=tmp_path / "ref.fa",
        reference_index=tmp_path / "ref.fa.fai",
        prefix="sample",
    )

    commands = stage.create_commands(Context(cpus=4))

    assert commands[1].processes[0].command == [
        "sort",
        "-k1,1",
        "-k2,2n",
        "sample.report.regions.bed",
    ]
    assert commands[1].processes[1].command == ["bedtools", "merge", "-i", "-"]
    assert commands[1].output_file == Path("sample.report.regions.merged.bed")
    assert commands[2].command == [
        "samtools",
        "index",
        str(tmp_path / "sample.bam"),
        "sample.report.alignment.index",
    ]
    assert commands[3].processes[0].command == [
        "samtools",
        "view",
        "--threads",
        "4",
        "-h",
        "-M",
        "-X",
        "-T",
        str(tmp_path / "ref.fa"),
        "-L",
        "sample.report.regions.merged.bed",
        str(tmp_path / "sample.bam"),
        "sample.report.alignment.index",
    ]
    assert commands[3].processes[1].command == [
        sys.executable,
        "-m",
        "snippy_ng",
        "utils",
        "aln",
        "samcrop",
        "--bed",
        "sample.report.regions.merged.bed",
    ]
    assert commands[3].processes[2].command == [
        "samtools",
        "sort",
        "--threads",
        "4",
        "-T",
        "sample.report.sort.tmp",
        "-O",
        "cram,level=9",
        "--reference",
        str(tmp_path / "ref.fa"),
        "-o",
        "sample.report.cram",
        "-",
    ]
    assert commands[4].command == [
        "samtools",
        "index",
        "sample.report.cram",
        "sample.report.cram.crai",
    ]


def test_sample_report_commands_without_alignment_only_render_table(tmp_path):
    stage = SampleReport(
        vcf=tmp_path / "sample.vcf",
        prefix="sample",
    )

    commands = stage.create_commands(Context(cpus=4))

    assert len(commands) == 2
    assert commands[0].description == "Create sample-report variant table and alignment windows"
    assert commands[1].description == "Render sample HTML report"
    assert commands[1].args[3] is None
    assert commands[1].args[4] is None
    assert commands[1].args[6] is None
    assert commands[1].args[7] is None
