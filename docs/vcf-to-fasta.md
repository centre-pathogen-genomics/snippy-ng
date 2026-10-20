---
title: One VCF file to rule them all
---

# One VCF file to rule them all

The `snippy.all.vcf.gz` produced by the `asm`, `short`, and `long` pipelines
carries everything needed to rebuild the masked consensus FASTA
(`snippy.fna`) with a few `bcftools`/`bedtools` commands. Beyond the real
variant calls, the all-VCF encodes depth regions as synthetic records:

| Kind | Meaning | Encoding |
| --- | --- | --- |
| `ZERODEPTH` | No aligned reads in the region | Explicit deletion: `REF=<bases>`, `ALT="-"`, `FILTER=PASS;ZERODEPTH`, `ID=DEL_n`, `INFO=ZERODEPTH;SVTYPE=DEL;END=...;SVLEN=-...` |
| `LOWDEPTH` span | Region below the min depth threshold (fixed `--min-depth` or auto `median − 2·MAD`) | N-substitution: `REF=<bases>`, `ALT="N"×len`, `FILTER=LOWDEPTHSPAN`, `ID=LOWDEPTH_n`, `INFO=LOWDEPTH;END=...` |

Example records:

```text
ctg	11	DEL_1	G	-	.	PASS;ZERODEPTH	TYPE=INDEL;ZERODEPTH;SVTYPE=DEL;END=11;SVLEN=-1	GT	1/1
ctg	21	LOWDEPTH_1	ACGTA	NNNNN	.	LOWDEPTHSPAN	LOWDEPTH;END=25	GT	1/1
```

The `LOWDEPTHSPAN` filter deliberately keeps the spans out of `snippy.pass.vcf`
(the pipeline's own consensus input) while `PASS;ZERODEPTH` keeps the
deletions in it. Both kinds of record are required to rebuild the FASTA, so
the recipe below selects them explicitly.

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
bcftools view -i '(FILTER~"PASS" || FILTER="LOWDEPTHSPAN") && (strlen(ALT)<=strlen(REF) || ALT="<DEL>")' "$VCF" \
  | bgzip -c > rec.vcf.gz
bcftools index -f rec.vcf.gz
bcftools consensus -f "$REF" --mark-del - --mark-snv lc -o rec.fa rec.vcf.gz
IN_FA=rec.fa   # current FASTA; each optional step below may update it
```

What each part of the filter does:

- `FILTER~"PASS"` — PASS variants *and* `PASS;ZERODEPTH` deletion blocks.
  `FILTER="PASS"` would be exact-match and silently drop the ZERODEPTH
  records.
- `FILTER="LOWDEPTHSPAN"` — the N-spans. They are intentionally not PASS, so
  `bcftools view -f PASS` would drop them.
- `strlen(ALT)<=strlen(REF) || ALT="<DEL>"` — drop insertions, whose
  breakpoint geometry is ambiguous for consensus. This is the same
  pre-consensus filter the pipeline applies.
- `--mark-del -` and `--mark-snv lc` — deleted bases become `-` and applied
  SNVs are written lowercase, matching the pipeline's pseudo-alignment.

Filtered variant calls are not applied. Low-depth calls that fall inside a
span are covered by the span's substitution; a call flagged `LOWDEPTH` whose
actual per-base position depth is not below the threshold sits outside any
span and is simply excluded by `FILTER="LOWDEPTHSPAN"`.

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

### 3. Uppercase the low-depth Ns (optional)

The consensus marks the spans lowercase (`n` via `--mark-snv lc`) but the
pipeline writes uppercase `N` at those positions:

```console
bcftools query -i 'INFO/LOWDEPTH=1' -f '%CHROM\t%POS0\t%INFO/END\n' "$VCF" \
  | sort -k1,1 -k2,2n | bedtools merge -i - > lowdepth.bed
bedtools maskfasta -fi "$IN_FA" -bed lowdepth.bed -fo snippy.reconstructed.fna -fullHeader -mc N
```

(`$IN_FA` is `rec.mixed.fa` after step 2, or `rec.fa` if you skipped it.)

## When the result matches the pipeline

This recipe replays the pipeline's `BcftoolsPseudoAlignment`,
`MaskMixedSites`, and `ApplyDepthMaskToFasta` stages on the all-VCF, so the
reconstructed FASTA is byte-identical to the pipeline output whenever the
reconstruction inputs agree with the pipeline's:

- Both the LOWDEPTH spans and the pipeline's depth mask derive from the same
  `DepthBedsFromBam` `min_depth_bed`, so their intervals should coincide; the
  reconstruction tests guard this invariant.
- ZERODEPTH deletions are `PASS;ZERODEPTH`, so they are applied both by the
  pipeline's pass-VCF consensus and by step 1.
- MIXEDSITE masking uses the same intervals and the same lowercase `n`.

Caveats:

- A user-supplied `--mask` BED (the `ApplyMask` step) is an external file and
  cannot be derived from the VCF. If you used `--mask`, apply that BED with
  `bedtools maskfasta -mc X` after step 3.
- Steps 2 and 3 only adjust masking case (`n` vs `N`). For the raw
  pseudo-alignment, step 1 alone is sufficient.