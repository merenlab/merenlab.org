---
layout: page
title: interactive [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DISPLAY.png" alt="DISPLAY" style="width:100px; border:none" />

A DISPLAY-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-display-contigs-stats](../../programs/anvi-display-contigs-stats)</span> <span class="artifact-p">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-p">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-p">[anvi-display-structure](../../programs/anvi-display-structure)</span> <span class="artifact-p">[anvi-inspect](../../programs/anvi-inspect)</span> <span class="artifact-p">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-p">[anvi-script-snvs-to-interactive](../../programs/anvi-script-snvs-to-interactive)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This page describes general properties of anvi'o interactive displays and programs that offer anvi'o interactive artifacts.

## Terminology

Anvi'o uses a simple terminology to address various aspects of interactive displays it produces, such as items, layers, views, orders, and so on. The purpose of this section is to provide some insights into these terminology using the figure below:

![an anvi'o display](../../images/anvio_display_template.png){:.center-img}

Even though the figure is a product of <span class="artifact-n">[anvi-display-pan](/software/anvio/help/7/programs/anvi-display-pan)</span>, the general terminology does not change across different interfaces, including the default visualizations of <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span>. Here are the descriptions of numbered areas in the figure:

* The tree denoted by **(1)** shows the organization of each `item`. Items could be contigs, gene clusters, bins, genes, or anything else depending on which mode the anvi'o interactive interface was initiated. The structure that orders items and denoted by **(1)** in the figure can be a phylogenetic or phylogenomic tree, or a dendrogram produced by a hierarchical clustering algorithm. In addition, there may be nothing there, if the user has requested or set a linear items order through <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7/artifacts/misc-data-items-order)</span>.
* Each concentric circle underneath the number **(2)** is called a `layer` and the data shown for items and layers as a whole is called a `view`. A **layer** can be a genome, a metagenome, or anything else depending on which mode the anvi'o interactive was initiated. The **view** is like a data table where a datum is set for each **item** in each **layer**. The view data is typically computed by anvi’o and stored in pan databases by <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/7/programs/anvi-pan-genome)</span> or profile databases by <span class="artifact-n">[anvi-profile](/software/anvio/help/7/programs/anvi-profile)</span>. The user add another view to the relevant combo box in the interface by providing a TAB-delimited file to <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span> through the command line argument `--additional-view`, or add new layers to extend these vies with additional data through <span class="artifact-n">[misc-data-items](/software/anvio/help/7/artifacts/misc-data-items)</span>.
* The tree denoted by **(3)** shows a specific ordering of layers. Anvi'o will compute various layer orders automatically based on available **view** depending on the analysis or visualization mode, and users can extend available **layer orders** through <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7/artifacts/misc-data-layer-orders)</span>.
* What is shown by **(4)** is the additional data for layers. the user can extend this section with additional information on layers using the <span class="artifact-n">[misc-data-layers](/software/anvio/help/7/artifacts/misc-data-layers)</span>.

The orchestrated use of <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7/programs/anvi-import-misc-data)</span>, <span class="artifact-n">[anvi-export-misc-data](/software/anvio/help/7/programs/anvi-export-misc-data)</span>, and <span class="artifact-n">[anvi-delete-misc-data](/software/anvio/help/7/programs/anvi-delete-misc-data)</span> provides a powerful framework to decorate items or layers in a display and enhance visualization of complex data. Please take a look at the following article on how to extend anvi'o displays:

* [https://merenlab.org/2017/12/11/additional-data-tables/](https://merenlab.org/2017/12/11/additional-data-tables/)

## Programs that give interactive access

If you're new to the anvi'o interactive interface, you'll probably want to check out [this tutorial for beginners](http://merenlab.org/tutorials/interactive-interface/) or the other resources on the  <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span> page. 

However, there are more interfaces available in anvi'o than just that one, so let's list them out: 

- <span class="artifact-n">[anvi-display-structure](/software/anvio/help/7/programs/anvi-display-structure)</span> lets you examine specific protein structures, along with SCV and SAAVs within it. (It even has [its own software page.](http://merenlab.org/software/anvio-structure/). It's kind of a big deal.)

- <span class="artifact-n">[anvi-display-contigs-stats](/software/anvio/help/7/programs/anvi-display-contigs-stats)</span> shows you various stats about the contigs within a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>, such as their hmm-hits, lengths, N and L statistics, and so on.

- <span class="artifact-n">[anvi-display-metabolism](/software/anvio/help/7/programs/anvi-display-metabolism)</span> is still under development but will allow you to interactively view metabolism estimation data using <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> under the hood. 

- <span class="artifact-n">[anvi-display-pan](/software/anvio/help/7/programs/anvi-display-pan)</span> displays information about the gene clusters that are stored in a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span>. It lets you easily view your core and accessory genes, and can even be turned into a metapangenome through importing additional data tables. 

- <span class="artifact-n">[anvi-inspect](/software/anvio/help/7/programs/anvi-inspect)</span> lets you look at a single split across your samples, as well as the genes identified within it. This interface can also be opened from the <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span> interface by asking for details about a specific split. 

- <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span> displays the information in a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. It lets you view the distribution of your contigs across your samples, manually bin metagenomic data into MAGSs (and refine those bins with <span class="artifact-n">[anvi-refine](/software/anvio/help/7/programs/anvi-refine)</span>), and much more. You can also use this to look at your genes instead of your contigs or [examine the genomes after a phylogenomic anlysis](http://merenlab.org/2017/06/07/phylogenomics/). Just look at that program page for a glimpse of this program's amazingness. 

- <span class="artifact-n">[anvi-script-snvs-to-interactive](/software/anvio/help/7/programs/anvi-script-snvs-to-interactive)</span> lets you view a comprehensive summary of the SNVs, SCVs, and SAAVs within your contigs. 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/interactive.md) to update this information.

