---
layout: page
title: genes-db [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/genes-db
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-gen-gene-level-stats-databases](../../programs/anvi-gen-gene-level-stats-databases)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-db-info](../../programs/anvi-db-info)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-migrate](../../programs/anvi-migrate)</span> <span class="artifact-r">[anvi-search-sequence-motifs](../../programs/anvi-search-sequence-motifs)</span></p>


## Description

An anvi'o genes database is a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>-like database that contains statistics, such their coverage and detection across samples, rather than contigs in a given <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>.

A gene database for a given <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> stored in a <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> will be automatically generated when <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> is run in 'gene mode'. For details, see the [relevant section](../programs/anvi-interactive/#visualizing-genes-instead-of-contigs) in <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span>

Alternatively, genes databases can be explicitly generated using the program <span class="artifact-n">[anvi-gen-gene-level-stats-databases](/software/anvio/help/7.1/programs/anvi-gen-gene-level-stats-databases)</span>. By default, this program will generate a gene database for each <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> for a given <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>. 

Due to the strucutral similarities between a <span class="artifact-n">[genes-db](/software/anvio/help/7.1/artifacts/genes-db)</span> and a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>, many of the anvi'o programs that operate on profile databases will also run on genes databases. These programs include those that import/export states and import/export misc additional data.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/genes-db.md) to update this information.

