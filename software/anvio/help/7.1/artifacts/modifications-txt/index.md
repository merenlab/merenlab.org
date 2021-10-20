---
layout: page
title: modifications-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/modifications-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-tabulate-trnaseq](../../programs/anvi-tabulate-trnaseq)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-plot-trnaseq](../../programs/anvi-plot-trnaseq)</span></p>


## Description

This tabular file contains data on predicted modifications in tRNA-seq seeds.

This file is produced by <span class="artifact-n">[anvi-tabulate-trnaseq](/software/anvio/help/7.1/programs/anvi-tabulate-trnaseq)</span>. The artifact for that program describes this and related tables in detail.

This tab-delimited file can be easily manipulated by the user. It is required input for <span class="artifact-n">[anvi-plot-trnaseq](/software/anvio/help/7.1/programs/anvi-plot-trnaseq)</span>.

## Example

The modifications shown in this table are from the seeds represented in the <span class="artifact-n">[seeds-specific-txt](/software/anvio/help/7.1/artifacts/seeds-specific-txt)</span> and <span class="artifact-n">[seeds-non-specific-txt](/software/anvio/help/7.1/artifacts/seeds-non-specific-txt)</span> example tables.

| gene_callers_id | contig_name | anticodon | aa | domain | phylum | class | order | family | genus | species | taxon_percent_id | seed_position | ordinal_name | ordinal_position | canonical_position | reference | sample_name | A | C | G | T |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | 19 | d_loop_beta_1 | 22 | 20 | G | DB_01 | 142 | 589 | 69411 | 1315 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | 19 | d_loop_beta_1 | 22 | 20 | G | DB_03 | 217 | 1056 | 83751 | 2592 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | 19 | d_loop_beta_1 | 22 | 20 | G | DB_05 | 42 | 212 | 28784 | 515 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | 19 | d_loop_beta_1 | 22 | 20 | G | DB_07 | 102 | 429 | 45633 | 977 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 32 | anticodon_loop_1 | 36 | 32 | T | DB_01 | 0 | 51 | 14 | 77 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 32 | anticodon_loop_1 | 36 | 32 | T | DB_03 | 1 | 274 | 97 | 642 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 32 | anticodon_loop_1 | 36 | 32 | T | DB_05 | 0 | 78 | 17 | 137 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 32 | anticodon_loop_1 | 36 | 32 | T | DB_07 | 0 | 19 | 18 | 87 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 37 | anticodon_loop_6 | 41 | 37 | G | DB_01 | 0 | 1 | 137 | 5 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 37 | anticodon_loop_6 | 41 | 37 | G | DB_03 | 5 | 18 | 916 | 64 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 37 | anticodon_loop_6 | 41 | 37 | G | DB_05 | 6 | 3 | 222 | 7 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | 37 | anticodon_loop_6 | 41 | 37 | G | DB_07 | 0 | 15 | 104 | 1 |


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/modifications-txt.md) to update this information.

