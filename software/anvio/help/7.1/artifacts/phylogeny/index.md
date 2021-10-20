---
layout: page
title: phylogeny [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/phylogeny
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/NEWICK.png" alt="NEWICK" style="width:100px; border:none" />

A NEWICK-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-items-order](../../programs/anvi-export-items-order)</span> <span class="artifact-p">[anvi-gen-phylogenomic-tree](../../programs/anvi-gen-phylogenomic-tree)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-items-order](../../programs/anvi-import-items-order)</span> <span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-script-checkm-tree-to-interactive](../../programs/anvi-script-checkm-tree-to-interactive)</span></p>


## Description

This is a NEWICK-formatted tree that describes the phylogenic relationships of your data. 

{:.notice}
Wondering what the NEWICK format is? Then you're in luck! It has its own [Wikipedia page](https://en.wikipedia.org/wiki/Newick_format).

### How to get one of these? 

You can use <span class="artifact-n">[anvi-gen-phylogenomic-tree](/software/anvio/help/7.1/programs/anvi-gen-phylogenomic-tree)</span> to create a phylogeny based on a series of genes. 

As discussed on the page for <span class="artifact-n">[anvi-gen-phylogenomic-tree](/software/anvio/help/7.1/programs/anvi-gen-phylogenomic-tree)</span>, you can also use an external program to get a NEWICK-formatted tree and use that. 

### What can I do with it? 

Firstly, you can use it to reorder elements of the interactive interface. To import this to rearrange the orders that your items appear (in other words, as the central phylogenetic tree when you open the interface), import it using <span class="artifact-n">[anvi-import-items-order](/software/anvio/help/7.1/programs/anvi-import-items-order)</span>. To import this as a tree describing your layers (the concentric circles in the anvi'o interface), convert this to a <span class="artifact-n">[misc-data-layer-orders-txt](/software/anvio/help/7.1/artifacts/misc-data-layer-orders-txt)</span> and use the program <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>.

Secondly, as done in the [Phylogenetics tutorial](http://merenlab.org/2017/06/07/phylogenomics/#working-with-fasta-files), you can open it in the interactive interface without an associated <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. To do this, run <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> as so:

<div class="codeblock" markdown="1">
anvi&#45;interactive &#45;t <span class="artifact&#45;n">[phylogeny](/software/anvio/help/7.1/artifacts/phylogeny)</span> \
                 &#45;&#45;title "Phylogenomics Tutorial" \
                 &#45;&#45;manual
</div>

This will create an empty <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> to store any <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>s you create and other such data. You can also add various information, such as taxonomy hits, as done in that same [Phylogenetics tutorial](http://merenlab.org/2017/06/07/phylogenomics/#working-with-fasta-files). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/phylogeny.md) to update this information.

