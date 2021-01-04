---
layout: page
title: misc-data-layer-orders [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-delete-misc-data](../../programs/anvi-delete-misc-data)</span> <span class="artifact-r">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span></p>


## Description

This is the section of your <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>/<span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> that contains custom additional information about the order that your layers are displayed in and the tree that relates them to each other . When you run <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span>, this data will determine what order the concentric circles are displayed in, as well as the tree that appears above the <span class="artifact-n">[misc-data-layers](/software/anvio/help/7/artifacts/misc-data-layers)</span> graphs.

As also defined in [this blog post](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology), this type of data will include information about how your layers are related to each other and determines their order. This data is stored as a tree that is displayed in the top-right. 

This data can be imported from a <span class="artifact-n">[misc-data-layer-orders-txt](/software/anvio/help/7/artifacts/misc-data-layer-orders-txt)</span> artifact with <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7/programs/anvi-import-misc-data)</span>.  It is also displayed when you run <span class="artifact-n">[anvi-show-misc-data](/software/anvio/help/7/programs/anvi-show-misc-data)</span> and can be exported or deleted with <span class="artifact-n">[anvi-export-misc-data](/software/anvio/help/7/programs/anvi-export-misc-data)</span> and <span class="artifact-n">[anvi-delete-misc-data](/software/anvio/help/7/programs/anvi-delete-misc-data)</span> respectively. 

For example, you could use this tree to indicate and group together samples that came from the same geographic location, samples that came from the same donor, samples of the same type,  samples collected with the same collection method, and so on. 

This is also used to import the taxonomy information at the end of [the pangenomics and phylogenomics workflow](http://merenlab.org/2017/06/07/phylogenomics/#pangenomic--phylogenomics). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-layer-orders.md) to update this information.

