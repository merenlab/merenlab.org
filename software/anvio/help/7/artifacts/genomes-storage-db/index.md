---
layout: page
title: genomes-storage-db [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-gen-genomes-storage](../../programs/anvi-gen-genomes-storage)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-analyze-synteny](../../programs/anvi-analyze-synteny)</span> <span class="artifact-r">[anvi-compute-functional-enrichment](../../programs/anvi-compute-functional-enrichment)</span> <span class="artifact-r">[anvi-compute-gene-cluster-homogeneity](../../programs/anvi-compute-gene-cluster-homogeneity)</span> <span class="artifact-r">[anvi-db-info](../../programs/anvi-db-info)</span> <span class="artifact-r">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-r">[anvi-get-sequences-for-gene-calls](../../programs/anvi-get-sequences-for-gene-calls)</span> <span class="artifact-r">[anvi-get-sequences-for-gene-clusters](../../programs/anvi-get-sequences-for-gene-clusters)</span> <span class="artifact-r">[anvi-meta-pan-genome](../../programs/anvi-meta-pan-genome)</span> <span class="artifact-r">[anvi-migrate](../../programs/anvi-migrate)</span> <span class="artifact-r">[anvi-pan-genome](../../programs/anvi-pan-genome)</span> <span class="artifact-r">[anvi-search-functions](../../programs/anvi-search-functions)</span> <span class="artifact-r">[anvi-split](../../programs/anvi-split)</span> <span class="artifact-r">[anvi-summarize](../../programs/anvi-summarize)</span> <span class="artifact-r">[anvi-update-db-description](../../programs/anvi-update-db-description)</span></p>


## Description

This is an Anvi'o database that **stores information about your genomes, primarily for use in pangenomic analyses.**

You can think of it like this: in a way, a genomes-storage-db is to the [the pangenomic workflow](http://merenlab.org/2016/11/08/pangenomics-v2/#generating-an-anvio-genomes-storage) what a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> is to the [the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). They both describe key information unique to your particular dataset and are required to run the vast majority of programs. 

### What kind of information? 

A genomes storage database contains information about the genomes that you inputted to create it, as well as the genes within them. 

Specifically, there are three tables stored in a genomes storage database: 

* A table describing the information about each of your genomes, such as their name, type (internal or external), GC content, number of contigs, completition, redunduncy, number of genes, etc. 
* A table describing the genes within your genomes. For each gene, this includes its gene caller id, associated genome and position, sequence, length, and whether or not it is partial. 
* A table describing the functions of your genes, including their sources and e-values. 

### Cool. How do I make one? 

You can generate one of these from an <span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span> (genomes described in <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span>s), <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span> (genomes described in <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>s), or both using the program <span class="artifact-n">[anvi-gen-genomes-storage](/software/anvio/help/7/programs/anvi-gen-genomes-storage)</span>. 

### Cool cool. What can I do with one? 

With one of these, you can run <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/7/programs/anvi-pan-genome)</span> to get a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span>. If a genomes storage database is the <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> of pangenomics, then a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> is the <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. It contains lots of information that is vital for analysis, and most programs will require both the <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> and its genomes storage database as an input. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/genomes-storage-db.md) to update this information.

