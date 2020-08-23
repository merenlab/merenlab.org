---
layout: post
title: misc-data-layers [artifact]
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

This is the section of your <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span>/<span class="artifact-n">[pan-db](/software/anvio/help/artifacts/pan-db)</span> that contains custom additional information about each of the layers of the interactive interface (usually displayed as the concentric circles). When you run <span class="artifact-n">[anvi-interactive](/software/anvio/help/programs/anvi-interactive)</span>, this data will appear as additional graphs in line with your layers, similar to how the sample names are displayed at the top. 

As also defined in [this blog post](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology), this type of data will include information about each layer of the interface (usually representing your samples). This data can be either numerical or categorical and can be imported into another database by the user through the command <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/programs/anvi-import-misc-data)</span>. 

If you would like to change the order that your layers are displayed, take a look at <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/artifacts/misc-data-layer-orders)</span>. Or, if you want to specifically import taxnomic information for your layers (if applicable), check out <span class="artifact-n">[anvi-import-taxonomy-for-layers](/software/anvio/help/programs/anvi-import-taxonomy-for-layers)</span>.

For example, this information could describe the salinity of a series of ocean samples, the continent your samples were taken in, or which of several collection methods was used. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-layers.md) to update this information.

