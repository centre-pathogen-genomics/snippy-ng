import json

from snippy_ng.stages.reporting import EpiReport


def test_epi_report_preserves_small_branch_lengths():
    report = EpiReport()
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
