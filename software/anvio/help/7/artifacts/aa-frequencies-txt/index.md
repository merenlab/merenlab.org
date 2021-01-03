---
layout: page
title: aa-frequencies-txt [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-get-aa-counts](../../programs/anvi-get-aa-counts)</span> <span class="artifact-p">[anvi-get-codon-frequencies](../../programs/anvi-get-codon-frequencies)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This file contains **the frequency of each amino acid for some reference context in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>**.  

This is a tab-delimited table where each column represents an amino acid and each row represents a specific reference context (most often this will be a gene after running <span class="artifact-n">[anvi-get-codon-frequencies](/software/anvio/help/7/programs/anvi-get-codon-frequencies)</span>). The numbers will either refer to counts of each amino acid or precent normalizations depending on the parameters with which you ran <span class="artifact-n">[anvi-get-codon-frequencies](/software/anvio/help/7/programs/anvi-get-codon-frequencies)</span>. 

You can also use <span class="artifact-n">[anvi-get-aa-counts](/software/anvio/help/7/programs/anvi-get-aa-counts)</span> to get this information for a <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span>, <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>, or <span class="artifact-n">[splits-txt](/software/anvio/help/7/artifacts/splits-txt)</span>. 

### Example

    gene_caller_id  Ala Arg Thr Asp ...
        1           0   0   1   2
        2           1   0   0   2
        .
        .
        .


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/aa-frequencies-txt.md) to update this information.

