---
layout: page
title: genes-db [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-gen-gene-level-stats-databases](../../programs/anvi-gen-gene-level-stats-databases)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-db-info](../../programs/anvi-db-info)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-migrate](../../programs/anvi-migrate)</span></p>


## Description

A genes database is a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> that also contains all of the gene-level statistics (such as coverage and detection). This is the output of the program <span class="artifact-n">[anvi-gen-gene-level-stats-databases](/software/anvio/help/7/programs/anvi-gen-gene-level-stats-databases)</span>, which will generate one of these per <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span> in the <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> you input. 

As such, you can use this as an input anywhere a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> is asked for. See that page for more information. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/genes-db.md) to update this information.

