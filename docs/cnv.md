# `snippy-ng utils cnv`

Estimate copy number variation from an aligned BAM or CRAM.

The command uses `samtools coverage` to calculate mean depth for each contig. By default, The largest contig is assumed to have one copy number. Other contigs are assigned rounded integer copy numbers from their depth relative to that baseline.


## Contig Copy Number

Estimate copy number for each contig:

```console
snippy-ng utils cnv sample.cram > sample.cnv.tsv
```

Output:

```text
contig_id	read_depth	copy_number
```

## Feature Copy Number

Supply a GFF to estimate copy number for each selected feature. The default
feature type is `gene`.

```console
snippy-ng utils cnv sample.cram --gff reference.gff > sample.features.cnv.tsv
```

Use another GFF feature type:

```console
snippy-ng utils cnv sample.cram --gff reference.gff --feature gene > sample.genes.cnv.tsv
```

Feature mode uses `samtools depth -aa -b` and reports the median per-base
depth for each feature.

Output:

```text
feature_id	contig_id	start	end	read_depth	copy_number
```

## Known Single-Copy Baseline

Use a known single-copy interval on the largest contig as the baseline:

```console
snippy-ng utils cnv sample.cram --known-single-copy 4518,5000 > sample.cnv.tsv
```

Use an explicit contig:

```console
snippy-ng utils cnv sample.cram --known-single-copy chr1:4518-5000 > sample.cnv.tsv
```

The known single-copy interval can also be used with feature mode:

```console
snippy-ng utils cnv sample.cram \
  --gff reference.gff \
  --known-single-copy 4518,5000 \
  > sample.features.cnv.tsv
```

The baseline is the median per-base depth across the known single-copy region.

## Options

| Option | Purpose |
| --- | --- |
| `ALIGNMENT` | Input BAM or CRAM. |
| `--gff PATH` | Estimate copy number per selected GFF feature instead of per contig. |
| `--feature TEXT` | GFF feature type to use with `--gff`; default `gene`. |
| `--known-single-copy TEXT` | Baseline interval as `START,END` on the largest contig or `CONTIG:START-END`. |
| `--no-header` | Do not print the TSV header. |
