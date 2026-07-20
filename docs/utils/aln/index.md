---
title: Alignment utilities
---

# Alignment utilities

Alignment utilities operate on SAM, BAM, or CRAM files produced by the mapping
stages.

## Commands

| Command | Purpose |
| --- | --- |
| [samclip](samclip.md) | Filter clipped reads from SAM alignments. |
| [cnv](cnv.md) | Estimate contig or feature copy numbers from aligned BAM or CRAM depth. |

## Notes

The samclip command is boundary-aware. A read is filtered when either terminal
clip exceeds the threshold, except when the clipped end reaches the contig
boundary. Fraction mode uses the total terminal clipping divided by the
original query length, including hard clips.