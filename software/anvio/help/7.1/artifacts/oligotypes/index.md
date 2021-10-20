---
layout: page
title: oligotypes [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/oligotypes
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-oligotype-linkmers](../../programs/anvi-oligotype-linkmers)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

Oligotyping is a computational strategy that **partitions a given set of sequences into homogeneous groups using only a subset of target nucleotide positions**. 

### History

Oligotyping was [first described in 2011](https://doi.org/10.1371/journal.pone.0026732) and has primarily been applied to 16S ribosomal RNA gene amplicons to resolve closely related but distinct that differ as low as one nucleotide at the amplified region, exceeding the sensitivity of the popular strategy of the time, 97% OTU clustering. Other papers that demonstrate the strengths of this approach include the following: [1](https://doi.org/10.1111/2041-210X.12114), [2](https://doi.org/10.1073/pnas.1409644111), [3](https://doi.org/10.1038/ismej.2014.195), and [4](https://doi.org/10.3389/fmicb.2014.00568).

The following figure from [the oligotyping methods paper](https://doi.org/10.1111/2041-210X.12114) depicts major steps of an oligotyping analysis for amplicon sequences:

![Oligotyping](../../images/oligotyping.jpg)

Although, Shannon entropy is not the only approach to identify highly variable nucleotide positions of interest, and they an also be provided by the user.

### Applications to metagenomics

This strategy also applies to metagenomic sequences that are mapped to a genomic context to describe the diversity of variable regions that are fully covered by short reads. In the context of metagenomic read recruitment, variable nucleotide positions can be chosen by the user from the positions of single-nucleotide variants [anvi'o recovers](https://merenlab.org/2015/07/20/analyzing-variability/) and presents through inspection pages in the interactive interface or through <span class="artifact-n">[variability-profile](/software/anvio/help/7.1/artifacts/variability-profile)</span>. An example application of oligotyping to metagenomics is demonstrated here:

* [An application of oligotyping in the metagenomic context: Oligotyping AmoC](https://merenlab.org/2015/12/09/musings-over-commamox/#an-application-of-oligotyping-in-the-metagenomic-context-oligotyping-amoc)

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/oligotypes.md) to update this information.

