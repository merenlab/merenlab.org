---
layout: page
title: pn-ps-data [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/pn-ps-data
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-get-pn-ps-ratio](../../programs/anvi-get-pn-ps-ratio)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This describes the output of <span class="artifact-n">[anvi-get-pn-ps-ratio](/software/anvio/help/7.1/programs/anvi-get-pn-ps-ratio)</span>, which calculates the pN/pS ratio for each gene in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. 

{:.notice}
See the page for <span class="artifact-n">[anvi-get-pn-ps-ratio](/software/anvio/help/7.1/programs/anvi-get-pn-ps-ratio)</span> for an explanation of the pN/pS ratio 

This describes a directory that contains the following four files: 

`pNpS.txt`: a long-format table of the pN/pS values, along with the groupby variables:

|   | corresponding_gene_call | sample_id   | pNpS_reference       |
| - | ----------------------- | ----------- | -------------------- |
| 0 | 1744                    | ANE_004_05M | 0.043503524536208836 |
| 1 | 1744                    | ANE_004_40M | 0.043628712253629943 |
| 2 | 1744                    | ANE_150_05M | 0.03810623760551494  |
| 3 | 1744                    | ANE_150_40M | 0.040815421982026576 |

`pN.txt`: a long-format table of the pN values, along with the groupby variables:

|   | corresponding_gene_call | sample_id   | pN_reference       |
| - | ----------------------- | ----------- | ------------------ |
| 0 | 1744                    | ANE_004_05M | 11.827627600424583 |
| 1 | 1744                    | ANE_004_40M | 11.106801744995472 |
| 2 | 1744                    | ANE_150_05M | 9.62355553228605   |
| 3 | 1744                    | ANE_150_40M | 10.067364489809782 |

`pS.txt`: a long-format table of the pS values, along with the groupby variables:

|   | corresponding_gene_call | sample_id   | pS_reference       |
| - | ----------------------- | ----------- | ------------------ |
| 0 | 1744                    | ANE_004_05M | 271.87745651689016 |
| 1 | 1744                    | ANE_004_40M | 254.57551165909962 |
| 2 | 1744                    | ANE_150_05M | 252.54541348089631 |
| 3 | 1744                    | ANE_150_40M | 246.6558962502711  |

`num_SCVs.txt`: a long-format table of the number of SCVs belonging to each group:

|   | corresponding_gene_call | sample_id   | num_SCVs |
| - | ----------------------- | ----------- | -------- |
| 0 | 1744                    | ANE_004_05M | 180      |
| 1 | 1744                    | ANE_004_40M | 166      |
| 2 | 1744                    | ANE_150_05M | 162      |
| 3 | 1744                    | ANE_150_40M | 160      |



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/pn-ps-data.md) to update this information.

