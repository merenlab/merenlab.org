---
layout: page
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


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This artifact describes the **order of items in visualization tasks**.

In anvi'o, main display items (such as 'gene clusters' in a pan database, 'contigs' in a profile database, etc) can be ordered either by a NEWICK formatted tree (such as a phylogenetic tree or a hierarchical clustering dendrogram), or by an array (such as a flat list of item names).

When a NEWICK tree is used to order items, it will appear as the tree in the central section of the default anvi'o interactive interface. When a flat list of items are provided to order items, the central display where a tree appears will be blank and the displayed items will still be ordered according to the list. In order words, items order is to <span class="artifact-n">[misc-data-items](/software/anvio/help/7/artifacts/misc-data-items)</span> as <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7/artifacts/misc-data-layer-orders)</span> is to <span class="artifact-n">[misc-data-layers](/software/anvio/help/7/artifacts/misc-data-layers)</span>: a description not of the items themselves, but of what order they go in on the interface. 

Anvi'o programs such as <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/7/programs/anvi-pan-genome)</span>, <span class="artifact-n">[anvi-merge](/software/anvio/help/7/programs/anvi-merge)</span>, and <span class="artifact-n">[anvi-profile](/software/anvio/help/7/programs/anvi-profile)</span> automatically generate NEWICK-formatted items order if possible (i.e., if you have less than 20,000 items). When you run these programs, they will put this information into your resulting <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> or <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. 

You can also export this information to give to a fellow anvi'o user or import this information if you have your own phylogenetic tree or desired order for your contigs.

You can use <span class="artifact-n">[anvi-import-items-order](/software/anvio/help/7/programs/anvi-import-items-order)</span> to import specific orders for your items, or <span class="artifact-n">[anvi-export-items-order](/software/anvio/help/7/programs/anvi-export-items-order)</span> to export this information.

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-items-order.md) to update this information.

