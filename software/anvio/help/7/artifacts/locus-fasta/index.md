---
layout: page
title: locus-fasta [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/FASTA.png" alt="FASTA" style="width:100px; border:none" />

A FASTA-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-locus](../../programs/anvi-export-locus)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

A locus-fasta is one of the outputs of <span class="artifact-n">[anvi-export-locus](/software/anvio/help/7/programs/anvi-export-locus)</span>, which creates exports specific regions of interest out of a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>. 

This artifact specifically describes the <span class="artifact-n">[fasta](/software/anvio/help/7/artifacts/fasta)</span> file that contains the sequence of one of the hits to the locus. 

This file is contained within the directory specified by the `-o` parameter and is named with the prefix defined by the `-O` parameter, followed by a numerical identifier for this particular hit. The sequence in this fasta file is also contained in the <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> of the same name. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/locus-fasta.md) to update this information.

