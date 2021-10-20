---
layout: page
title: trna-taxonomy [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/trna-taxonomy
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-trna-taxonomy](../../programs/anvi-run-trna-taxonomy)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-estimate-trna-taxonomy](../../programs/anvi-estimate-trna-taxonomy)</span></p>


## Description

This contains the taxonomic annotations for each of the tRNA sequences found in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, which are the results of running <span class="artifact-n">[anvi-run-trna-taxonomy](/software/anvio/help/7.1/programs/anvi-run-trna-taxonomy)</span>. 

You can use this information to estimate the taxnomy of genomes, metagenomes, or collections stored in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> using the program <span class="artifact-n">[anvi-estimate-trna-taxonomy](/software/anvio/help/7.1/programs/anvi-estimate-trna-taxonomy)</span>

Recall that this information was calculated using [GTDB](https://gtdb.ecogenomic.org/), so it might not be entirely accurate for Eukaryotic tRNAs. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/trna-taxonomy.md) to update this information.

