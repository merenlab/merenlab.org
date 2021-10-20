---
layout: page
title: genes-stats [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/genes-stats
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/STATS.png" alt="STATS" style="width:100px; border:none" />

A STATS-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-gen_stats_for_single_copy_genes.py](../../programs/anvi-script-gen_stats_for_single_copy_genes.py)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This file contains information about your genes. 

It is a tab-delimited text file where each row represents a specific gene and each column provides different information. 

As of now, the only program that returns data in this format is <span class="artifact-n">[anvi-script-gen_stats_for_single_copy_genes.py](/software/anvio/help/7.1/programs/anvi-script-gen_stats_for_single_copy_genes.py)</span>, which returns this information for the single copy core genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. 

From left to right, these tell you 
* The source for this gene (ex `Protista_83`)
* The name of the contig that this gene is a part of
* The gene name 
* The e-value (of the HMM hit that was used to find this gene)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/genes-stats.md) to update this information.

