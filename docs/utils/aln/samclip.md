---
title: Samclip
---

# `snippy-ng utils aln samclip`

Filter clipped reads from a SAM file using a boundary-aware terminal clipping
rule.

The command checks the terminal clipping at both ends of each alignment and
filters the read when either end exceeds the configured threshold. If a clipped
end reaches the contig boundary, that end is ignored. This means reads with
large clipped segments can remain when the clipping is explained by the contig
edge rather than an internal alignment artifact.

`--max-clip-fraction` uses the total terminal clipping divided by the original
query length, including hard-clipped bases.

## Example

Filter clipped alignments with the default maximum terminal clip length:

```console
snippy-ng utils aln samclip --index ref.fa.fai input.sam > clipped.sam
```

Use only the fraction-based rule:

```console
snippy-ng utils aln samclip \
  --index ref.fa.fai \
  --max-clip-fraction 0.1 \
  input.sam > clipped.sam
```

Combine the absolute and fraction-based rules:

```console
snippy-ng utils aln samclip \
  --index ref.fa.fai \
  --max 10 \
  --max-clip-fraction 0.1 \
  input.sam > clipped.sam
```

## Options

| Option | Purpose |
| --- | --- |
| `SAM_FILE` | Input SAM file. Reads from standard input when omitted. |
| `--index PATH` | Reference FASTA index (.fai) for the reference used during alignment. |
| `--fix-mate / --no-fix-mate` | Attempt to mark mates as unmapped when one read in a pair is filtered. |
| `--max INT` | Maximum terminal clip length to allow on either end, after ignoring clipping at contig boundaries. |
| `--max-clip-fraction FLOAT` | Maximum total terminal clipping as a fraction of the original read length, including hard clips. |
| `--invert` | Keep only clipped reads instead of filtering them out. |
| `--debug` | Write debug information about clipped reads to standard error. |

## Behavior notes

- `--max` is applied to each end independently.
- A read is filtered when the left end, the right end, or both exceed the
  threshold.
- Clipping at the contig boundary is ignored before the comparison.
- `--max-clip-fraction` measures total terminal clipping, not just soft clips.
- `--fix-mate` expects queryname-sorted SAM input.