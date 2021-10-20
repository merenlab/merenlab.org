---
layout: page
title: hmm-hits [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/hmm-hits
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-hmms](../../programs/anvi-run-hmms)</span> <span class="artifact-p">[anvi-scan-trnas](../../programs/anvi-scan-trnas)</span> <span class="artifact-p">[anvi-script-filter-hmm-hits-table](../../programs/anvi-script-filter-hmm-hits-table)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-delete-hmms](../../programs/anvi-delete-hmms)</span> <span class="artifact-r">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span> <span class="artifact-r">[anvi-run-scg-taxonomy](../../programs/anvi-run-scg-taxonomy)</span> <span class="artifact-r">[anvi-script-filter-hmm-hits-table](../../programs/anvi-script-filter-hmm-hits-table)</span> <span class="artifact-r">[anvi-script-gen-hmm-hits-matrix-across-genomes](../../programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span> <span class="artifact-r">[anvi-script-get-hmm-hits-per-gene-call](../../programs/anvi-script-get-hmm-hits-per-gene-call)</span></p>


## Description

The search results for an <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span> in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. Essentially, this is the part of a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> that handles the HMM data. In anvi'o, this is usually functional annotations, such as identifying specfic ribosomal RNAs, various single-copy core genes, and transfer RNAs, though the user can also define their own HMM sources. 

Upon creation, a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> will not contain any HMM results. In order to populate it, users can run <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7.1/programs/anvi-run-hmms)</span> using any <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>. The program <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/7.1/programs/anvi-scan-trnas)</span> also populates a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>'s hmm-hits with potential tranfer RNA hits.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/hmm-hits.md) to update this information.

