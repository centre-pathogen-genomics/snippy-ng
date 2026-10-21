---
title: One VCF file to rule them all
---

# One VCF file to rule them all

The `snippy.all.vcf.gz` produced by the `asm`, `short`, and `long` pipelines
carries everything needed to rebuild the masked consensus FASTA
(`snippy.fna`) with a few `bcftools`/`bedtools` commands. Beyond the real
variant calls, the all-VCF encodes depth regions as synthetic records:

| Kind | Meaning | Encoding | Rebuilt as |
| --- | --- | --- | --- |
| `ZERODEPTH` | No aligned reads in the region | Left-anchored symbolic deletion block: `POS=<base before the region>`, `REF=<that base>`, `ALT="<DEL>"`, `FILTER=ZERODEPTH`, `ID=DEL_n`, `INFO=TYPE=INDEL;ZERODEPTH;SVTYPE=DEL;END=...;SVLEN=-...` | `-` |
| `LOWDEPTH` span | Region below the min depth threshold (fixed `--min-depth` or auto `median − 2·MAD`) | The same symbolic deletion block, marked with `FILTER=LOWDEPTHSPAN`, `ID=LOWDEPTH_n`, `INFO=TYPE=LOWDEPTH;LOWDEPTH;SVTYPE=DEL;END=...;SVLEN=-...` | `N` |

Example records:

```text
Wildtype	35187	DEL_2	C	<DEL>	.	ZERODEPTH	TYPE=INDEL;ZERODEPTH;SVTYPE=DEL;END=35747;SVLEN=-560	GT	1/1
Wildtype	20	LOWDEPTH_1	G	<DEL>	.	LOWDEPTHSPAN	TYPE=LOWDEPTH;LOWDEPTH;SVTYPE=DEL;END=35187;SVLEN=-35167	GT	1/1
```

Both kinds are `<DEL>` blocks that the consensus step applies as deletions, so
the depth regions come out as `-`; step 3 then masks the low-depth intervals
with `N`, which is what the pipeline's own depth mask does. Because the block
is symbolic, `REF` holds a single anchor base instead of the whole region —
for the two records above that is the difference between ~35 kb of `REF`/`ALT`
text and a handful of bytes.

A symbolic block needs a base before the region to anchor on, so two cases
still use an explicit `REF` sequence:

- a region that starts at position 1 of a contig — `POS=1`,
  `REF=<region bases>`, `ALT="-"`, `END=<last base of the region>`;
- a region that covers a whole contig — `POS=1`, `REF=<whole contig>`,
  `ALT="-"`, `END=<contig length>`. This is not strict VCF, but it preserves
  deletion semantics for consensus generation.

Neither kind of block is `PASS`, so both are excluded from `snippy.pass.vcf`
and only appear in the all-VCF. Both are required to rebuild the FASTA, so the
recipe below selects them explicitly. The pipeline's own consensus runs on a
filtered view of the all-VCF (`PASS` records plus the `ZERODEPTH` blocks), so
zero-depth regions still become `-`; it gets its `N`s from the depth BED,
which step 3 reproduces from the VCF.

`asm` runs differ: a default `asm` run (no input `--vcf`, `add_deletions_to_vcf`
enabled) adds ZERODEPTH deletions covering regions of the reference with no
assembled sequence, but no LOWDEPTH spans, because assembly-based calling has
no read-depth threshold. Only steps 1 and 2 apply to `asm` output; step 3 is a
no-op there. If `asm` was given an input `--vcf` or deletion synthesis was
disabled, its all-VCF has no ZERODEPTH blocks either.

## Rebuild the FASTA

With the output directory as the working directory:

```console
REF=reference/ref.fa     # the prepared reference; the name matches your input reference
VCF=snippy.all.vcf.gz
```

### 1. Consensus

```console
bcftools view -i '(FILTER~"PASS" || FILTER="ZERODEPTH" || FILTER="LOWDEPTHSPAN") && (strlen(ALT)<=strlen(REF) || ALT="<DEL>")' "$VCF" \
  | bgzip -c > rec.vcf.gz
bcftools index -f rec.vcf.gz
bcftools consensus -f "$REF" --mark-del - --mark-snv lc -o rec.fa rec.vcf.gz
IN_FA=rec.fa   # current FASTA; each optional step below may update it
```

What each part of the filter does:

- `FILTER~"PASS"` — the real PASS variant calls.
- `FILTER="ZERODEPTH"` — the zero-depth blocks, which are deliberately not
  PASS so that they stay out of `snippy.pass.vcf`.
- `FILTER="LOWDEPTHSPAN"` — the low-depth blocks. They are intentionally not
  PASS either, so `bcftools view -f PASS` would drop both block kinds.
- `strlen(ALT)<=strlen(REF) || ALT="<DEL>"` — drop insertions, whose
  breakpoint geometry is ambiguous for consensus, while keeping the symbolic
  `<DEL>` depth blocks (their `ALT` is longer than their one-base `REF`). This
  is the same pre-consensus filter the pipeline applies.
- `--mark-del -` and `--mark-snv lc` — deleted bases become `-` and applied
  SNVs are written lowercase, matching the pipeline's pseudo-alignment. Both
  depth block kinds are deletions here, so both become `-`; step 3 turns the
  low-depth runs into `N`.

Filtered variant calls are not applied. Low-depth calls that fall inside a
span are covered by the span's deletion (and then masked in step 3); a call
flagged `LOWDEPTH` whose actual per-base position depth is not below the
threshold sits outside any span and is simply excluded by
`FILTER="LOWDEPTHSPAN"`.

### 2. Mask MIXEDSITE sites (optional)

Flagged `MIXEDSITE` records are real calls, not depth regions, so they are
excluded from the consensus and masked instead, exactly as the pipeline does:

```console
bcftools query -i 'FILTER="MIXEDSITE"' -f '%CHROM\t%POS0\t%REF\n' "$VCF" \
  | awk '{print $1"\t"$2"\t"($2+length($3))}' \
  | sort -k1,1 -k2,2n | bedtools merge -i - > mixed.bed
bedtools maskfasta -fi "$IN_FA" -bed mixed.bed -fo rec.mixed.fa -fullHeader -mc n
IN_FA=rec.mixed.fa
```

The intervals are sorted before merging exactly as the pipeline's
`MaskMixedSites` stage does, so the recipe mirrors the pipeline
command-for-command.

### 3. Mask the low-depth regions with Ns

The consensus applies the low-depth blocks as deletions, but the pipeline
writes `N` at those positions, so mask the same intervals again:

```console
bcftools query -i 'INFO/LOWDEPTH=1' -f '%CHROM\t%POS\t%INFO/END\n' "$VCF" \
  | sort -k1,1 -k2,2n | bedtools merge -i - > lowdepth.bed
bedtools maskfasta -fi "$IN_FA" -bed lowdepth.bed -fo snippy.reconstructed.fna -fullHeader -mc N
```

(`$IN_FA` is `rec.mixed.fa` after step 2, or `rec.fa` if you skipped it.)

`%POS` is the block's anchor base, one before the region it removes, which is
exactly the 0-based start of the low-depth interval. Using `%POS0` would start
the mask one base early and overwrite the anchor base.

## When the result matches the pipeline

This recipe replays the pipeline's `BcftoolsPseudoAlignment`,
`MaskMixedSites`, and `ApplyDepthMaskToFasta` stages on the all-VCF, so the
reconstructed FASTA is byte-identical to the pipeline output whenever the
reconstruction inputs agree with the pipeline's:

- Both the LOWDEPTH spans and the pipeline's depth mask derive from the same
  `DepthBedsFromBam` `min_depth_bed`, so their intervals should coincide: the
  interval `%POS`–`%INFO/END` from a span is the same 0-based interval the BED
  holds. The reconstruction tests guard this invariant.
- ZERODEPTH blocks are not PASS, so the pipeline's consensus runs on a
  filtered view of the all-VCF (PASS records plus those blocks); step 1
  selects the same records.
- MIXEDSITE masking uses the same intervals and the same lowercase `n`.

Caveats:

- A user-supplied `--mask` BED (the `ApplyMask` step) is an external file and
  cannot be derived from the VCF. If you used `--mask`, apply that BED with
  `bedtools maskfasta -mc X` after step 3.
- Steps 2 and 3 are both needed for the masked consensus: without step 2
  MIXEDSITE bases stay unmasked, and without step 3 low-depth bases stay `-`
  instead of `N`. For the raw, unmasked pseudo-alignment, step 1 alone is
  enough.