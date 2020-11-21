---
layout: page
title: structure-db [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-gen-structure-database](../../programs/anvi-gen-structure-database)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-3dev](../../programs/anvi-3dev)</span> <span class="artifact-r">[anvi-db-info](../../programs/anvi-db-info)</span> <span class="artifact-r">[anvi-export-structures](../../programs/anvi-export-structures)</span> <span class="artifact-r">[anvi-gen-fixation-index-matrix](../../programs/anvi-gen-fixation-index-matrix)</span> <span class="artifact-r">[anvi-gen-variability-profile](../../programs/anvi-gen-variability-profile)</span> <span class="artifact-r">[anvi-update-structure-database](../../programs/anvi-update-structure-database)</span></p>

## Description

This ddatabase contains the **predicted 3D structure data** of the sequences encoded by contigs in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. It is the result of running <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/programs/anvi-gen-structure-database)</span> (which uses MODELLER to predict your protein strucutres initially based on alignment to sequences with known strucutres). 

You can use this database to generate a variability profile (with <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/programs/anvi-gen-variability-profile)</span>) and then use that data (along with this database) to visualize your 3D structures and their relation to SCVs and SAAVs with <span class="artifact-n">[anvi-3dev](/software/anvio/help/programs/anvi-3dev)</span>. 

Besides this, you can export the structure data into an external `.pdb` file (using <span class="artifact-n">[anvi-export-structures](/software/anvio/help/programs/anvi-export-structures)</span>) or generation the fixation index matrix (with <span class="artifact-n">[anvi-gen-fixation-index-matrix](/software/anvio/help/programs/anvi-gen-fixation-index-matrix)</span>). 

For more information on the structure database, see [this blog post](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#the-structure-database). 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/structure-db.md) to update this information.

