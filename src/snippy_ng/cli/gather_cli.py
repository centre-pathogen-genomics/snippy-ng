import click

from snippy_ng.cli.utils import absolute_path_callback


@click.command(context_settings={'show_default': True})
@click.argument("inputs", required=False, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, nargs=-1)
@click.option("--reference", "--ref", required=False, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Reference genome to include in JSON output and exclude from the search")
@click.option("--max-depth", "-d", type=click.INT, default=4, help="Maximum directory depth to search for sequence files", show_default=True)
@click.option("--exclude", "-e", type=click.STRING, default=None, help="Regular expression to exclude files based on their name", show_default=True)
@click.option("--aggressive-ids", "-a", is_flag=True, default=False, help="Aggressively parse sample IDs from file paths", show_default=True)
@click.option("--json", "-j", is_flag=True, default=False, help="Output JSON instead of TSV", show_default=True)
def gather(**config):
    """
    Gather samples into a CSV file for multi-sample analysis

    Examples:

        $ snippy-ng utils gather > samples.csv
    """
    from snippy_ng.utils.gather import gather_samples_config
    import os

    gathered = gather_samples_config(
        inputs=config["inputs"] if config.get("inputs") else [os.getcwd()],
        max_depth=config["max_depth"],
        aggressive_ids=config["aggressive_ids"],
        exclude_name_regex=config["exclude"],
        reference=config.get("reference"),
    )

    if config["json"]:
        import json
        print(json.dumps(gathered, indent=2))
    else:
        import csv
        import sys
        # Collect all possible inner keys
        samples = gathered["samples"]
        fieldnames = set()
        for sample_data in samples.values():
            fieldnames.update(sample_data.keys())

        # Deterministic order, with sample first
        fieldnames.discard("type")
        fieldnames = ["sample", "type"] + sorted(fieldnames)

        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()

        for sample, sample_data in samples.items():
            row = {"sample": sample}
            row.update(sample_data)
            writer.writerow(row)
        

