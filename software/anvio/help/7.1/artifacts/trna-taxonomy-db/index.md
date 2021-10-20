---
layout: page
title: trna-taxonomy-db [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/trna-taxonomy-db
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-setup-trna-taxonomy](../../programs/anvi-setup-trna-taxonomy)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-run-trna-taxonomy](../../programs/anvi-run-trna-taxonomy)</span></p>


## Description

This artifact represents the [GTDB](https://gtdb.ecogenomic.org/) data (from [Parks et al. 2018](https://doi.org/10.1038/nbt.4229)) downloaded by <span class="artifact-n">[anvi-setup-trna-taxonomy](/software/anvio/help/7.1/programs/anvi-setup-trna-taxonomy)</span>. This information is required to run <span class="artifact-n">[anvi-run-trna-taxonomy](/software/anvio/help/7.1/programs/anvi-run-trna-taxonomy)</span> and <span class="artifact-n">[anvi-estimate-trna-taxonomy](/software/anvio/help/7.1/programs/anvi-estimate-trna-taxonomy)</span>. 

{:.notice}
If the results from this tRNA taxonomy search end up in a paper, make sure to cite [Parks et al. 2018](https://doi.org/10.1038/nbt.4229) for their information.

By default, it is stored at `anvio/data/misc/TRNA-TAXONOMY`. This directory contains a few files for each anticodon, each forming a fancy search database so that you can associate tRNA reads in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> with taxonomy information. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/trna-taxonomy-db.md) to update this information.

