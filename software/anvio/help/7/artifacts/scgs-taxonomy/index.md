---
layout: page
title: scgs-taxonomy [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-scg-taxonomy](../../programs/anvi-run-scg-taxonomy)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-estimate-scg-taxonomy](../../programs/anvi-estimate-scg-taxonomy)</span></p>


## Description

This contains the taxonomy annotions for each of the single-copy core genes found in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>; in other words, this contains the results of a run of <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/7/programs/anvi-run-scg-taxonomy)</span>. 

This information was found through the [GTDB](https://gtdb.ecogenomic.org/) database, so it will not work with Eukaryotic genomes. 

This information allows you to quickly estimate the taxonomy of genomes, metagenomes, or bins stored in your contigs-db by using the command <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/scgs-taxonomy.md) to update this information.

