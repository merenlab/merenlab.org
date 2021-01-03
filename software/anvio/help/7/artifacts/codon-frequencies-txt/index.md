---
layout: page
title: codon-frequencies-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-get-codon-frequencies](../../programs/anvi-get-codon-frequencies)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This file contains **the frequency of each codon in each gene in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>.** 

This is a tab-delimited table where each column represents a codon and each row represents a specific gene. The numbers will either refer to counts of each codon or precent normalizations depending on the parameters with which you ran <span class="artifact-n">[anvi-get-codon-frequencies](/software/anvio/help/7/programs/anvi-get-codon-frequencies)</span>. 

### Example

    gene_caller_id  GCA GCC GCG GCT ...
        1           0   0   1   2
        2           1   0   0   2
        .
        .
        .


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/codon-frequencies-txt.md) to update this information.

