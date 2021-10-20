---
layout: page
title: dendrogram [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/dendrogram
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/NEWICK.png" alt="NEWICK" style="width:100px; border:none" />

A NEWICK-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-experimental-organization](../../programs/anvi-experimental-organization)</span> <span class="artifact-p">[anvi-export-items-order](../../programs/anvi-export-items-order)</span> <span class="artifact-p">[anvi-matrix-to-newick](../../programs/anvi-matrix-to-newick)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-items-order](../../programs/anvi-import-items-order)</span> <span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span></p>


## Description

This described a [NEWICK-formatted](https://en.wikipedia.org/wiki/Newick_format) tree that is not representative of the phylogenic relationships between your samples. 

{:.notice}
If you're looking for phylogenic trees, take a look at <span class="artifact-n">[phylogeny](/software/anvio/help/7.1/artifacts/phylogeny)</span> 

Instead, the dendrogram artifact most often describes the tree used as a <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7.1/artifacts/misc-data-items-order)</span>: the order that the items in <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> are displayed in (the central tree in the circular display). Often, these are the order of your contigs or genes based on their relatedness to each other (for example by tetranucleotide frequency or differencial coverage). These trees are also contained in <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7.1/artifacts/misc-data-layer-orders)</span>.

A dendrogram is also listed as the output of programs that are not necessarily related to phylogenetics (like <span class="artifact-n">[anvi-matrix-to-newick](/software/anvio/help/7.1/programs/anvi-matrix-to-newick)</span>).  


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/dendrogram.md) to update this information.

