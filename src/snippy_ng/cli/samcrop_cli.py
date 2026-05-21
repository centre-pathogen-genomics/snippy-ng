import sys

import click

from snippy_ng.cli.utils import AbsolutePath


@click.command(context_settings={'show_default': True})
@click.argument("sam_file", required=False, nargs=1, type=AbsolutePath(exists=True, readable=True))
@click.option("--bed", "-b", required=True, type=AbsolutePath(exists=True, readable=True), help="BED file of intervals to retain")
def samcrop(sam_file, bed):
    """
    Hard-crop SAM alignments to intervals in a BED file.

    Examples:

        $ samtools view -h input.bam | snippy-ng utils aln samcrop --bed windows.bed | samtools view -b -o cropped.bam -
    """
    from snippy_ng.utils.samcrop import samcrop_filter_lines, load_bed

    intervals = load_bed(bed)
    with open(sam_file, "r", encoding="utf-8") if sam_file else sys.stdin as sam_lines:
        output_buffer = []
        for line in samcrop_filter_lines(sam_lines, intervals):
            output_buffer.append(line)
            if len(output_buffer) >= 8192:
                sys.stdout.write("".join(output_buffer))
                output_buffer.clear()
        if output_buffer:
            sys.stdout.write("".join(output_buffer))
