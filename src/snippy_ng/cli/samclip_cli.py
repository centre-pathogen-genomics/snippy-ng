import click

from snippy_ng.cli.utils import absolute_path_callback


@click.command(context_settings={'show_default': True})
@click.argument("sam_file", required=False, nargs=1, type=click.Path(exists=True, readable=True), callback=absolute_path_callback)
@click.option("--index", "-i", required=True, type=click.Path(exists=True, readable=True), callback=absolute_path_callback, help="Reference FASTA index (.fai) file corresponding to the reference used for alignment")
@click.option("--fix-mate/--no-fix-mate", default=True, help="Attempt to fix mate information for paired-end reads when one read is filtered out", show_default=True)
@click.option("--max", "-m", type=click.INT, default=10, help="Maximum clip length to allow before filtering out the read", show_default=True)
@click.option("--invert", is_flag=True, default=False, help="Invert the filter to keep only clipped reads instead of filtering them out")
@click.option("--debug", is_flag=True, default=False, help="Output debug information about clipped reads to stderr")
def samclip(
    sam_file,
    index,
    max,
    fix_mate,
    invert,
    debug,
):
    """
    Filter clipped reads from a SAM file

    Examples:

        $ snippy-ng utils samclip --index ref.fa.fai input.sam > clipped.sam
    """
    import sys
    from snippy_ng.utils.samclip import samclip_filter_lines, fai_to_dict

    # Load reference index
    with open(index, 'r') as f:
        contig_lengths = fai_to_dict(f)
    # Run the pipeline
    with open(sam_file, 'r') if sam_file else sys.stdin as sam_lines:
        for line in samclip_filter_lines(
            sam_lines,
            contig_lengths=contig_lengths,
            max_clip=max,
            invert=invert,
            on_debug=lambda msg: click.echo(msg, err=True) if debug else None,
            fix_mate=fix_mate,
        ):
            print(line, end="")

    


    