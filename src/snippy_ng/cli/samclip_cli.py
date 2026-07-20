import click

from snippy_ng.cli.utils import AbsolutePath


@click.command(context_settings={'show_default': True})
@click.argument("sam_file", required=False, nargs=1, type=AbsolutePath(exists=True, readable=True))
@click.option("--index", "-i", required=True, type=AbsolutePath(exists=True, readable=True), help="Reference FASTA index (.fai) file corresponding to the reference used for alignment")
@click.option("--fix-mate/--no-fix-mate", default=True, help="Attempt to fix mate information for paired-end reads when one read is filtered out", show_default=True)
@click.option("--max", "-m", type=click.INT, default=10, help="Maximum terminal clip length to allow on either end before filtering, after ignoring clips at contig boundaries", show_default=True)
@click.option(
    "--max-clip-fraction",
    type=click.FloatRange(min=0, max=1),
    default=None,
    help="Maximum total terminal clipping as a fraction of the original read length, including hard clips and after ignoring boundary clips",
)
@click.option("--invert", is_flag=True, default=False, help="Invert the filter to keep only clipped reads instead of filtering them out")
@click.option("--debug", is_flag=True, default=False, help="Output debug information about clipped reads to stderr")
@click.pass_context
def samclip(
    ctx,
    sam_file,
    index,
    max,
    max_clip_fraction,
    fix_mate,
    invert,
    debug,
):
    """
    Filter clipped reads from a SAM file using boundary-aware terminal clipping.

    ``--max`` filters reads when either terminal clip exceeds the threshold,
    unless the clipped end reaches the contig boundary. ``--max-clip-fraction``
    filters reads by total terminal clipping fraction, using the original query
    length including hard clips.

    Examples:

        $ snippy-ng utils aln samclip --index ref.fa.fai input.sam > clipped.sam
    """
    import sys
    from snippy_ng.utils.samclip import samclip_filter_lines, fai_to_dict

    # Keep the historical default absolute threshold for existing invocations.
    # When a fraction is supplied on its own, it replaces that default; an
    # explicitly supplied --max combines with the fraction criterion.
    max_clip = max if (
        max_clip_fraction is None
        or ctx.get_parameter_source("max") != click.core.ParameterSource.DEFAULT
    ) else None

    # Load reference index
    with open(index, 'r') as f:
        contig_lengths = fai_to_dict(f)
    # Run the pipeline
    with open(sam_file, 'r') if sam_file else sys.stdin as sam_lines:
        output_buffer = []
        line_count = 0
        for line in samclip_filter_lines(
            sam_lines,
            contig_lengths=contig_lengths,
            max_clip=max_clip,
            max_clip_fraction=max_clip_fraction,
            invert=invert,
            on_debug=lambda msg: click.echo(message=msg, err=True) if debug else None,
            fix_mate=fix_mate,
        ):
            output_buffer.append(line)
            if len(output_buffer) >= 8192:
                sys.stdout.write("".join(output_buffer))
                output_buffer.clear()
            line_count += 1
            if line_count % 1000000 == 0:
                click.echo(f"[samclip] Processed {line_count} lines...", err=True)
        if output_buffer:
            sys.stdout.write("".join(output_buffer))
        click.echo(f"[samclip] Finished processing {line_count} lines.", err=True)
