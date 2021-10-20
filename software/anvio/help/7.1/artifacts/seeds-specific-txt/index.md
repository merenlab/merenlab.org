---
layout: page
title: seeds-specific-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/seeds-specific-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This tabular file contains data on the specific coverages of tRNA-seq seeds.

Specific coverage represents reads that are assigned uniquely to a tRNA seed. See the <span class="artifact-n">[trnaseq-profile-db](/software/anvio/help/7.1/artifacts/trnaseq-profile-db)</span> artifact for a fuller explanation of specific versus nonspecific coverage. The rows and columns of this table are identical to <span class="artifact-n">[seeds-non-specific-txt](/software/anvio/help/7.1/artifacts/seeds-non-specific-txt)</span> except the type of coverage data reported in each.

This file is produced by <span class="artifact-n">[anvi-tabulate-trnaseq](/software/anvio/help/7.1/programs/anvi-tabulate-trnaseq)</span>. The artifact for that program describes this and related tables in detail.

This tab-delimited file can be easily manipulated by the user. It is required input for <span class="artifact-n">[anvi-plot-trnaseq](/software/anvio/help/7.1/programs/anvi-plot-trnaseq)</span>.

## Example

The seeds shown in this table are also shown in the <span class="artifact-n">[seeds-non-specific-txt](/software/anvio/help/7.1/artifacts/seeds-non-specific-txt)</span> example. Modifications from these seeds are shown in the <span class="artifact-n">[modifications-txt](/software/anvio/help/7.1/artifacts/modifications-txt)</span> example.

| gene_callers_id | contig_name | anticodon | aa | domain | phylum | class | order | family | genus | species | taxon_percent_id | sample_name | mean_coverage | relative_mean_coverage | relative_discriminator_coverage | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 17a | 18 | 19 | 20 | 20a | 20b | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39 | 40 | 41 | 42 | 43 | 44.01 | 44.02 | 44.03 | 44.04 | 44.05 | 44.06 | 44.07 | 44.08 | 44.09 | 44.1 | 44.11 | 44.12 | 44.13 | 44.14 | 44.15 | 44.16 | 44.17 | 44.18 | 44.19 | 44.2 | 44.21 | 44.22 | 44.23 | 49 | 50 | 51 | 52 | 53 | 54 | 55 | 56 | 57 | 58 | 59 | 60 | 61 | 62 | 63 | 64 | 65 | 66 | 67 | 68 | 69 | 70 | 71 | 72 | 73 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_01 | 71456.2 | 0.25805309 | 0.2578848 |  | 71392 | 71398 | 71398 | 71400 | 71413 | 71414 | 71414 | 71414 | 71419 | 71425 | 71426 | 71426 | 71426 | 71437 | 71445 | 71451 | 71451 |  | 71451 | 71455 | 71457 |  |  | 71460 | 71467 | 71467 | 71471 | 71468 | 71470 | 71470 | 71470 | 71475 | 71489 | 71540 | 71534 | 71529 | 71606 | 71579 | 71582 | 71583 | 71586 | 71586 | 71583 | 71583 | 71583 | 71585 | 71584 | 71584 | 71584 | 71587 | 71587 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 71586 | 71586 | 71572 | 71572 | 71570 | 71504 | 71505 | 71504 | 71503 | 71503 | 71503 | 71503 | 71503 | 71503 | 71503 | 71503 | 71502 | 71500 | 71500 | 71500 | 71497 | 71496 | 71488 | 71466 | 68324 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_03 | 87746.6 | 0.22523292 | 0.22512637 |  | 87568 | 87582 | 87588 | 87594 | 87598 | 87599 | 87596 | 87597 | 87598 | 87606 | 87608 | 87611 | 87611 | 87611 | 87612 | 87615 | 87615 |  | 87614 | 87614 | 87616 |  |  | 87617 | 87619 | 87619 | 87620 | 87620 | 87620 | 87621 | 87619 | 87620 | 87630 | 87689 | 87683 | 87694 | 87732 | 87898 | 87912 | 87921 | 87925 | 87926 | 87926 | 87926 | 87931 | 87924 | 87925 | 87929 | 87981 | 87986 | 87986 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 87985 | 87985 | 87981 | 87980 | 87978 | 87952 | 87952 | 87951 | 87952 | 87954 | 87952 | 87950 | 87946 | 87946 | 87946 | 87943 | 87941 | 87939 | 87939 | 87939 | 87936 | 87935 | 87929 | 87909 | 84530 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_05 | 29533 | 0.22692849 | 0.22190602 |  | 29490 | 29494 | 29499 | 29509 | 29516 | 29516 | 29518 | 29521 | 29525 | 29528 | 29536 | 29536 | 29536 | 29538 | 29539 | 29547 | 29550 |  | 29550 | 29549 | 29553 |  |  | 29552 | 29552 | 29553 | 29552 | 29553 | 29553 | 29549 | 29550 | 29552 | 29553 | 29564 | 29561 | 29559 | 29638 | 29605 | 29602 | 29604 | 29607 | 29607 | 29607 | 29608 | 29607 | 29607 | 29607 | 29607 | 29607 | 29609 | 29609 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 29609 | 29609 | 29606 | 29606 | 29605 | 29601 | 29601 | 29601 | 29601 | 29601 | 29599 | 29599 | 29599 | 29599 | 29596 | 29596 | 29597 | 29595 | 29595 | 29594 | 29593 | 29504 | 29486 | 29451 | 26980 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_07 | 47065 | 0.181019 | 0.18087983 |  | 47078 | 47087 | 47105 | 47113 | 47114 | 47114 | 47114 | 47116 | 47124 | 47127 | 47129 | 47129 | 47131 | 47133 | 47139 | 47139 | 47139 |  | 47135 | 47138 | 47141 |  |  | 47145 | 47142 | 47147 | 47145 | 47145 | 47142 | 47143 | 47135 | 47134 | 47133 | 47154 | 47125 | 47093 | 47111 | 47051 | 47058 | 47067 | 47071 | 47073 | 47072 | 47069 | 47069 | 47069 | 47072 | 47072 | 47073 | 47073 | 47074 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 47074 | 47074 | 47073 | 47069 | 47069 | 47042 | 47043 | 47042 | 47038 | 47042 | 47042 | 47042 | 47042 | 47041 | 47041 | 47040 | 47040 | 47040 | 47040 | 47040 | 47033 | 47032 | 47031 | 47018 | 45352 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_01 | 142.4 | 0.00051437 | 0.00052465 |  | 141 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 |  | 142 | 142 | 142 | 142 |  | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 142 | 141 | 141 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 143 | 139 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_03 | 1006.4 | 0.00258316 | 0.00249815 |  | 1010 | 1010 | 1010 | 1010 | 1011 | 1011 | 1011 | 1011 | 1011 | 1010 | 1010 | 1010 | 1010 | 1010 | 1010 | 1010 | 1011 |  | 1011 | 1011 | 1014 | 1014 |  | 1014 | 1014 | 1014 | 1014 | 1014 | 1014 | 1014 | 1014 | 1014 | 1014 | 1014 | 1014 | 1013 | 1012 | 1012 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 1003 | 998 | 998 | 938 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_05 | 297.5 | 0.00228586 | 0.00285402 |  | 225 | 227 | 227 | 227 | 227 | 227 | 225 | 228 | 228 | 230 | 230 | 230 | 230 | 230 | 230 | 230 | 230 |  | 230 | 230 | 230 | 230 |  | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 232 | 233 | 238 | 360 | 367 | 367 | 367 | 370 | 370 | 370 | 370 | 370 | 370 | 370 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 370 | 362 | 362 | 347 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_07 | 120 | 0.00046169 | 0.00046664 |  | 116 | 116 | 116 | 116 | 116 | 116 | 116 | 116 | 116 | 119 | 119 | 119 | 118 | 118 | 118 | 118 | 118 |  | 119 | 122 | 126 | 126 |  | 126 | 126 | 126 | 126 | 126 | 126 | 126 | 126 | 126 | 126 | 126 | 124 | 124 | 122 | 123 | 120 | 120 | 120 | 120 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 119 | 117 | 117 | 117 |


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/seeds-specific-txt.md) to update this information.

