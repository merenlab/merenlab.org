---
layout: page
title: gene-taxonomy [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-import-taxonomy-for-genes](../../programs/anvi-import-taxonomy-for-genes)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This describes **the taxonomy information for the genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>**. 

You can use <span class="artifact-n">[anvi-import-taxonomy-for-genes](/software/anvio/help/7/programs/anvi-import-taxonomy-for-genes)</span> to import this information through a <span class="artifact-n">[gene-taxonomy-txt](/software/anvio/help/7/artifacts/gene-taxonomy-txt)</span>, either from external data or by using software like [Kaiju](https://github.com/bioinformatics-centre/kaiju) or [Centrifuge](https://github.com/infphilo/centrifuge). See [this blog post](http://merenlab.org/2016/06/18/importing-taxonomy/) for a comprehensive tutorial. 

Once this information is populated, it will be displayed in most downstream interfaces, including <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span>. 

You can also add taxonomy information for the layers in your interface (most likely sections of your sample when analyzing a single sample) using <span class="artifact-n">[anvi-import-taxonomy-for-layers](/software/anvio/help/7/programs/anvi-import-taxonomy-for-layers)</span>, or at the genome or metagenome level using <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/gene-taxonomy.md) to update this information.

