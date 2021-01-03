---
layout: page
title: bin [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/BIN.png" alt="BIN" style="width:100px; border:none" />

A BIN-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-cluster-contigs](../../programs/anvi-cluster-contigs)</span> <span class="artifact-p">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-p">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-p">[anvi-refine](../../programs/anvi-refine)</span> <span class="artifact-p">[anvi-rename-bins](../../programs/anvi-rename-bins)</span> <span class="artifact-p">[anvi-script-add-default-collection](../../programs/anvi-script-add-default-collection)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-estimate-scg-taxonomy](../../programs/anvi-estimate-scg-taxonomy)</span> <span class="artifact-r">[anvi-estimate-trna-taxonomy](../../programs/anvi-estimate-trna-taxonomy)</span> <span class="artifact-r">[anvi-gen-fixation-index-matrix](../../programs/anvi-gen-fixation-index-matrix)</span> <span class="artifact-r">[anvi-gen-gene-level-stats-databases](../../programs/anvi-gen-gene-level-stats-databases)</span> <span class="artifact-r">[anvi-gen-variability-profile](../../programs/anvi-gen-variability-profile)</span> <span class="artifact-r">[anvi-get-short-reads-from-bam](../../programs/anvi-get-short-reads-from-bam)</span> <span class="artifact-r">[anvi-get-split-coverages](../../programs/anvi-get-split-coverages)</span> <span class="artifact-r">[anvi-inspect](../../programs/anvi-inspect)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-merge-bins](../../programs/anvi-merge-bins)</span> <span class="artifact-r">[anvi-refine](../../programs/anvi-refine)</span> <span class="artifact-r">[anvi-rename-bins](../../programs/anvi-rename-bins)</span> <span class="artifact-r">[anvi-script-gen-distribution-of-genes-in-a-bin](../../programs/anvi-script-gen-distribution-of-genes-in-a-bin)</span></p>


## Description

A bin is, in its simplest form, **a group of contigs**.  (Think of a literal bin that you're putting data into.)

In Anvi'o, you'll most commonly work with bins both in this form and in the form of <span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span>, espeically when you want to work with bins contained in more than one <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. A group of bins is called a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>.

## What can you use bins for?

### Bins are all over 'omics
In general, you can also use bins and <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>s to limit what you are analyzing downstream. A ton of anvi'o programs are able to take in a bin or a collection (a group of bins) so that you don't have to analyze your entire <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> when you just want to look at one section. This is true of many, many anvi'o analyses you can run, from <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> to <span class="artifact-n">[anvi-get-sequences-for-hmm-hits](/software/anvio/help/7/programs/anvi-get-sequences-for-hmm-hits)</span>.

Below are some more specific examples of binning in action!

### Metagenomic binning
A common use of binning is in **metagenomics**. See [this tutorial](http://merenlab.org/2016/06/22/anvio-tutorial-v2/) for details, but essentially, the data you're working with in metagenomics is a sample of all of the genetic material in an environmental sample. It's like pulling a bunch of random DNA fragments out of a bucket of ocean water. If you want to try to rebuild the individual genomes from that mess, one common strategy is to piece together genomes de novo by trying to group the contigs together. This is called genome-resolved metagenomic binning.

Basically, in metagenomic binning, you're trying to group together a bunch of contigs that all belong to the same genome using various metrics like tetranucleotide frequency, differential coverage, completion, etc. You can do this either using various algorithms (for instance, those used by <span class="artifact-n">[anvi-cluster-contigs](/software/anvio/help/7/programs/anvi-cluster-contigs)</span>) or manually through the interactive interface (<span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span>). In this type of binning, a lot of your bins will be complete metagenome-assembled genomes, or MAGs; however, if you find an interesting group of contigs (for example, a prophage or a plasmid or even a particular domain), you can also put that into a bin. Then, you can group these bins in different ways using <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>s.

For more information, you might want to watch [this lovely lecture on genome-resolved metagenomics](https://www.youtube.com/watch?v=RjNdHGK4ruo).

### Pangenomic Workflows
You can also use bins to group together gene clusters. This is useful if you want a specific group of contigs to remain together through your entire analysis. Just provide your <span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span> file to <span class="artifact-n">[anvi-gen-genomes-storage](/software/anvio/help/7/programs/anvi-gen-genomes-storage)</span>.

Wow, this binning thing seems BINcredible! (not sorry)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/bin.md) to update this information.

