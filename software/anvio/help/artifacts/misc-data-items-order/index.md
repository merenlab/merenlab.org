---
layout: post
title: misc-data-items-order [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-import-items-order](../../programs/anvi-import-items-order)</span> <span class="artifact-p">[anvi-merge](../../programs/anvi-merge)</span> <span class="artifact-p">[anvi-pan-genome](../../programs/anvi-pan-genome)</span> <span class="artifact-p">[anvi-profile](../../programs/anvi-profile)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"></p>

## Description

This artifact describes the **order of the branches for a visualization**.

This is the tree in the central section of the default anvi'o interactive interface, or the order that your contigs (or genes or bins, depending on the mode) display in when in the circles. In order words, it is to <span class="artifact-n">[misc-data-items](/software/anvio/help/artifacts/misc-data-items)</span> as <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/artifacts/misc-data-layer-orders)</span> is to <span class="artifact-n">[misc-data-layers](/software/anvio/help/artifacts/misc-data-layers)</span>: a description not of the items themselves, but of what order they go in on the interface. 

Most often, this is a phylogenetic tree, but it can also just describe what order to display the branches or items in. 

As of now, this is an provided by programs that generate a tree of this kind, including <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/programs/anvi-pan-genome)</span>, <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span>, and <span class="artifact-n">[anvi-profile](/software/anvio/help/programs/anvi-profile)</span>. When you run these programs, they will put this information into your resulting <span class="artifact-n">[pan-db](/software/anvio/help/artifacts/pan-db)</span> or <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span>. 

However, you can also export this information to give to a fellow Anvi'o user or import this information if you have your own phylogenetic tree or desired order for your contigs. This is especially useful if you want to perform manual binning. 

To import your own order for your items, use <span class="artifact-n">[anvi-import-items-order](/software/anvio/help/programs/anvi-import-items-order)</span>. To export this information, use <span class="artifact-n">[anvi-export-items-order](/software/anvio/help/programs/anvi-export-items-order)</span>. 




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-items-order.md) to update this information.

