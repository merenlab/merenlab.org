---
layout: page
title: misc-data-items [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/misc-data-items
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-get-sequences-for-gene-clusters](../../programs/anvi-get-sequences-for-gene-clusters)</span> <span class="artifact-p">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-p">[anvi-search-sequence-motifs](../../programs/anvi-search-sequence-motifs)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-delete-misc-data](../../programs/anvi-delete-misc-data)</span> <span class="artifact-r">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span></p>


## Description

This is the section of your <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>/<span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> that contains custom additional information about each of the items in the central section of the interactive interface. When you run <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span>, this data will appear as additional concentric circles. 

As also defined in [this blog post](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology), this type of data will include information about each item (whether that's a contig, gene, or bin). This data is either numerical or categorical and can be imported into another database from a <span class="artifact-n">[misc-data-items-txt](/software/anvio/help/7.1/artifacts/misc-data-items-txt)</span> using <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>. It is also displayed when you run <span class="artifact-n">[anvi-show-misc-data](/software/anvio/help/7.1/programs/anvi-show-misc-data)</span> and can be exported or deleted with <span class="artifact-n">[anvi-export-misc-data](/software/anvio/help/7.1/programs/anvi-export-misc-data)</span> and <span class="artifact-n">[anvi-delete-misc-data](/software/anvio/help/7.1/programs/anvi-delete-misc-data)</span> respectively. 

To change the order that the items are displayed in, take a look at <span class="artifact-n">[anvi-import-items-order](/software/anvio/help/7.1/programs/anvi-import-items-order)</span>.

For example, this information could describe whether or not each bin reached a certain completion threshold, the e-score of the function annotation on each gene, or different categories that the total length of a contig could fall into (1-1.5 kb, 1.5-2 kb, 2-2.5 kb, and so on). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-items.md) to update this information.

