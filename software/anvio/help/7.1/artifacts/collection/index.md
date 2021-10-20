---
layout: page
title: collection [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/collection
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/COLLECTION.png" alt="COLLECTION" style="width:100px; border:none" />

A COLLECTION-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-cluster-contigs](../../programs/anvi-cluster-contigs)</span> <span class="artifact-p">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-p">[anvi-import-collection](../../programs/anvi-import-collection)</span> <span class="artifact-p">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-p">[anvi-rename-bins](../../programs/anvi-rename-bins)</span> <span class="artifact-p">[anvi-script-add-default-collection](../../programs/anvi-script-add-default-collection)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-cluster-contigs](../../programs/anvi-cluster-contigs)</span> <span class="artifact-r">[anvi-delete-collection](../../programs/anvi-delete-collection)</span> <span class="artifact-r">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-r">[anvi-estimate-genome-completeness](../../programs/anvi-estimate-genome-completeness)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-estimate-scg-taxonomy](../../programs/anvi-estimate-scg-taxonomy)</span> <span class="artifact-r">[anvi-estimate-trna-taxonomy](../../programs/anvi-estimate-trna-taxonomy)</span> <span class="artifact-r">[anvi-export-collection](../../programs/anvi-export-collection)</span> <span class="artifact-r">[anvi-gen-gene-level-stats-databases](../../programs/anvi-gen-gene-level-stats-databases)</span> <span class="artifact-r">[anvi-get-aa-counts](../../programs/anvi-get-aa-counts)</span> <span class="artifact-r">[anvi-get-split-coverages](../../programs/anvi-get-split-coverages)</span> <span class="artifact-r">[anvi-merge-bins](../../programs/anvi-merge-bins)</span> <span class="artifact-r">[anvi-rename-bins](../../programs/anvi-rename-bins)</span> <span class="artifact-r">[anvi-split](../../programs/anvi-split)</span> <span class="artifact-r">[anvi-summarize](../../programs/anvi-summarize)</span> <span class="artifact-r">[anvi-script-gen-distribution-of-genes-in-a-bin](../../programs/anvi-script-gen-distribution-of-genes-in-a-bin)</span> <span class="artifact-r">[anvi-script-gen-genomes-file](../../programs/anvi-script-gen-genomes-file)</span></p>


## Description

Essentially, a collection **is a group of <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>s**.

You can use collections to represent all kinds of things. The default collection (that you'll get by running <span class="artifact-n">[anvi-script-add-default-collection](/software/anvio/help/7.1/programs/anvi-script-add-default-collection)</span>) is called DEFAULT and contains all of your contigs. However, you can put any group of bins into their own collection, and use that to limit what you're analyzing downstream. A ton of anvi'o programs are able to take in a bin or a collection so that you don't have to analyze your entire contigs-db when you just want to look at one section. 

To look at the collections contained within an Anvi'o database, just run <span class="artifact-n">[anvi-show-collections-and-bins](/software/anvio/help/7.1/programs/anvi-show-collections-and-bins)</span>. Or, to look at the completion estimates for all bins within a collection, use <span class="artifact-n">[anvi-estimate-genome-completeness](/software/anvio/help/7.1/programs/anvi-estimate-genome-completeness)</span>.

To view the content within a collection, use <span class="artifact-n">[anvi-summarize](/software/anvio/help/7.1/programs/anvi-summarize)</span>.

### Examples of collections in action

If your bins represent MAGs, you could use a collection to group related MAGs together. For example, you could group together genomes from the same population, or from the same taxonomic order. Or your could go even wider and have a collection for all of the Archaea genomes in your sample. You're not limited by taxonomy either. You could go wild and have a collection for all of your bins that you suspect are prophages. 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/collection.md) to update this information.

