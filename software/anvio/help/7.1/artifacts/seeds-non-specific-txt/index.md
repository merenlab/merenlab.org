---
layout: page
title: seeds-non-specific-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/seeds-non-specific-txt
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

This tabular file contains data on the nonspecific coverages of tRNA-seq seeds.

Nonspecific coverage represents reads that are not unique to a single tRNA seed. See the <span class="artifact-n">[trnaseq-profile-db](/software/anvio/help/7.1/artifacts/trnaseq-profile-db)</span> artifact for a fuller explanation of specific versus nonspecific coverage. The rows and columns of this table are identical to <span class="artifact-n">[seeds-specific-txt](/software/anvio/help/7.1/artifacts/seeds-specific-txt)</span> except for the type of coverage data reported in each.

This file is produced by <span class="artifact-n">[anvi-tabulate-trnaseq](/software/anvio/help/7.1/programs/anvi-tabulate-trnaseq)</span>. The artifact for that program describes this and related tables in detail.

This tab-delimited file can be easily manipulated by the user. It is optional input for <span class="artifact-n">[anvi-plot-trnaseq](/software/anvio/help/7.1/programs/anvi-plot-trnaseq)</span>.

## Example

The seeds shown in this table are also shown in the <span class="artifact-n">[seeds-specific-txt](/software/anvio/help/7.1/artifacts/seeds-specific-txt)</span> example. Modifications from these seeds are shown in the <span class="artifact-n">[modifications-txt](/software/anvio/help/7.1/artifacts/modifications-txt)</span> example.

| gene_callers_id | contig_name | anticodon | aa | domain | phylum | class | order | family | genus | species | taxon_percent_id | sample_name | mean_coverage | relative_mean_coverage | relative_discriminator_coverage | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 17a | 18 | 19 | 20 | 20a | 20b | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39 | 40 | 41 | 42 | 43 | 44.01 | 44.02 | 44.03 | 44.04 | 44.05 | 44.06 | 44.07 | 44.08 | 44.09 | 44.1 | 44.11 | 44.12 | 44.13 | 44.14 | 44.15 | 44.16 | 44.17 | 44.18 | 44.19 | 44.2 | 44.21 | 44.22 | 44.23 | 49 | 50 | 51 | 52 | 53 | 54 | 55 | 56 | 57 | 58 | 59 | 60 | 61 | 62 | 63 | 64 | 65 | 66 | 67 | 68 | 69 | 70 | 71 | 72 | 73 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_01 | 100372.4 |  | 14497 | 14608 | 14815 | 14882 | 14985 | 15828 | 15854 | 15895 | 16410 | 16565 | 16840 | 16960 | 16975 | 16990 | 17490 | 18529 | 19087 |  | 19683 | 21763 | 24353 |  |  | 24699 | 25182 | 29097 | 30476 | 30609 | 30612 | 30491 | 30125 | 29973 | 29506 | 29417 | 26259 | 31169 | 145828 | 153750 | 155936 | 156187 | 156518 | 157233 | 157226 | 157178 | 157429 | 158124 | 159941 | 164453 | 167924 | 170595 | 170567 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 170572 | 170577 | 170547 | 170541 | 170326 | 169581 | 169509 | 169509 | 169497 | 169497 | 169494 | 169491 | 168727 | 168721 | 168719 | 168719 | 168719 | 168719 | 168719 | 168719 | 168489 | 168485 | 168480 | 167688 | 155628 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_03 | 203816.6 |  | 10498 | 10599 | 11105 | 11217 | 11255 | 11270 | 11350 | 11444 | 12349 | 12539 | 13028 | 13331 | 13337 | 13390 | 14325 | 15079 | 15603 |  | 15769 | 18168 | 21167 |  |  | 23927 | 24910 | 27041 | 28271 | 28395 | 28604 | 28612 | 28749 | 29242 | 30775 | 32254 | 33570 | 44299 | 319895 | 335328 | 337653 | 341382 | 342776 | 344543 | 345794 | 345762 | 345808 | 346281 | 347647 | 354052 | 360294 | 361948 | 361948 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 362317 | 362317 | 362303 | 362303 | 362226 | 362056 | 362056 | 362056 | 362056 | 362050 | 362046 | 362044 | 362044 | 362044 | 362044 | 362044 | 362032 | 362027 | 362027 | 362027 | 361349 | 361349 | 361349 | 360094 | 345770 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_05 | 26137.9 |  | 5111 | 5184 | 5259 | 5704 | 5979 | 5979 | 5979 | 6011 | 6550 | 6587 | 6587 | 6611 | 6611 | 6611 | 6936 | 7087 | 7090 |  | 7170 | 8243 | 9158 |  |  | 9488 | 9868 | 12268 | 12323 | 12866 | 12616 | 12506 | 12640 | 12630 | 12838 | 12292 | 11336 | 11621 | 36476 | 37479 | 38030 | 38892 | 39018 | 39018 | 39084 | 39068 | 39272 | 39272 | 39272 | 40666 | 41573 | 41879 | 41873 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 42019 | 42019 | 42015 | 42015 | 41996 | 41873 | 41857 | 41695 | 41689 | 41689 | 41689 | 41689 | 41683 | 41683 | 41683 | 41569 | 41569 | 40839 | 40839 | 40839 | 40839 | 40538 | 40495 | 40464 | 36174 |
0 | c_000000684460_DB_R05_06 | TAC | Val | Bacteria | Firmicutes | Clostridia | Lachnospirales | Lachnospiraceae |  |  | 100 | DB_07 | 182536.6 |  | 16048 | 16134 | 16358 | 16639 | 16664 | 16679 | 16679 | 16757 | 17351 | 17494 | 17494 | 17547 | 17613 | 17737 | 18346 | 18771 | 18776 |  | 19172 | 20831 | 21986 |  |  | 22549 | 22933 | 25352 | 26246 | 26370 | 25534 | 25798 | 25385 | 25277 | 25529 | 25758 | 25079 | 33202 | 289207 | 300224 | 303731 | 306231 | 306774 | 307692 | 307695 | 307604 | 307604 | 307828 | 308578 | 314195 | 317240 | 321017 | 321023 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 321634 | 321634 | 321615 | 321615 | 321596 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321066 | 321059 | 321059 | 320379 | 320373 | 320363 | 320230 | 303028 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_01 | 3923.1 |  | 32 | 32 | 32 | 32 | 32 | 32 | 32 | 32 | 32 | 32 | 33 | 33 | 33 | 33 | 33 | 33 | 33 |  | 33 | 50 | 79 | 82 |  | 82 | 82 | 82 | 82 | 82 | 82 | 82 | 90 | 90 | 90 | 96 | 109 | 118 | 117 | 125 | 163 | 320 | 7564 | 7948 | 7970 | 7991 | 7991 | 7991 | 8014 | 8014 | 8016 | 8041 | 8041 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 8041 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8046 | 8022 | 8022 | 8022 | 8021 | 7961 | 7955 | 7213 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_03 | 8502 |  | 27 | 27 | 27 | 29 | 29 | 31 | 31 | 43 | 49 | 59 | 59 | 59 | 59 | 59 | 59 | 59 | 59 |  | 59 | 60 | 65 | 65 |  | 65 | 65 | 65 | 65 | 65 | 65 | 65 | 72 | 72 | 78 | 78 | 92 | 116 | 140 | 173 | 236 | 679 | 16238 | 17144 | 17344 | 17428 | 17428 | 17452 | 17482 | 17482 | 17482 | 17482 | 17482 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 17482 | 17482 | 17482 | 17482 | 17482 | 17482 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17480 | 17411 | 17320 | 16198 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_05 | 1254.6 |  | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  | 0 | 4 | 18 | 18 |  | 18 | 18 | 18 | 18 | 26 | 26 | 26 | 32 | 32 | 32 | 32 | 40 | 45 | 45 | 60 | 62 | 89 | 2379 | 2492 | 2502 | 2562 | 2594 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2604 | 2557 | 2557 | 2557 | 2557 | 2446 | 2426 | 2058 |
1 | c_000000805276_DB_R05_05 | ACG | Arg | Bacteria | Firmicutes |  |  |  |  |  | 98.649 | DB_07 | 4217.8 |  | 60 | 60 | 60 | 60 | 60 | 60 | 60 | 64 | 72 | 78 | 78 | 78 | 78 | 78 | 78 | 78 | 78 |  | 78 | 90 | 107 | 107 |  | 107 | 113 | 113 | 113 | 113 | 119 | 119 | 117 | 117 | 117 | 120 | 120 | 132 | 140 | 172 | 212 | 410 | 7980 | 8475 | 8539 | 8553 | 8584 | 8594 | 8594 | 8594 | 8594 | 8599 | 8599 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8599 | 8546 | 8546 | 8124 |


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/seeds-non-specific-txt.md) to update this information.

