---
layout: page
title: split-bins [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/split-bins
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-split](../../programs/anvi-split)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This is the result of <span class="artifact-n">[anvi-split](/software/anvio/help/7.1/programs/anvi-split)</span>: self-contained anvi'o projects that contain just the contents of a single <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> from your original database. 

This describes a directory that either contains either a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> and <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> (if that's what you gave <span class="artifact-n">[anvi-split](/software/anvio/help/7.1/programs/anvi-split)</span> as an input) or a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> pair. The contigs or genomes and gene clusters described in these databases will be only those contained in the bin that the directory's name corresponds to.  


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/split-bins.md) to update this information.

