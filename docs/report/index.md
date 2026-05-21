---
title: Report
---

# Interactive reports

Snippy-NG provides two HTML report commands under `snippy-ng utils report`:

- [`tree`](tree.md): an interactive phylogenetic tree report from a Newick tree
- [`sample`](sample.md): an interactive per-sample variant report from a VCF, with optional alignment windows

## Commands

| Command | Purpose |
| --- | --- |
| `snippy-ng utils report tree` | Render an interactive phylogenetic tree report. |
| `snippy-ng utils report sample` | Render an interactive sample variant report. |

## Typical workflow

Generate a tree report after a multi-sample run:

```console
snippy-ng utils gather --json --ref tests/data/reference.gbk tests/data/{wildtype,mutant}* > samples.json
snippy-ng multi samples.json --cpus 6 -o multi
snippy-ng tree --fast multi/core/core.095.aln -o multi/tree
snippy-ng utils report tree multi/tree/tree.treefile --metadata multi/snippy.vcf.summary.tsv -o multi/report
```

Generate a sample report from an existing VCF:

```console
snippy-ng utils report sample short/snippy.pass.vcf -o short/report
```

See the dedicated pages for examples, screenshots, and command references.
