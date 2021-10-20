---
layout: page
title: hmm-hits-matrix-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/hmm-hits-matrix-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-gen-hmm-hits-matrix-across-genomes](../../programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This file is the output of <span class="artifact-n">[anvi-script-gen-hmm-hits-matrix-across-genomes](/software/anvio/help/7.1/programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span> and describes the <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span> across multiple genomes or bins for a single <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>. 

The first column describes each of the genomes (if the input was an <span class="artifact-n">[external-genomes](/software/anvio/help/7.1/artifacts/external-genomes)</span>) or bins (if the input was an <span class="artifact-n">[internal-genomes](/software/anvio/help/7.1/artifacts/internal-genomes)</span>) that the matrix describes. The following columns describe each of the genes in your <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>. The data within the table describes the number of hits that gene had in that genome or bin. 

For example, if you were to run <span class="artifact-n">[anvi-script-gen-hmm-hits-matrix-across-genomes](/software/anvio/help/7.1/programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span> with the `Bacteria_71` <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span> on two hypothetical genomes, you would get a file like this:

    genome_or_bin    ADK    AICARFT_IMPCHas    ATP-synt    ATP-synt_A    Adenylsucc_synt    Chorismate_synt    EF_TS    ...
    Genome_1         11     10                 9           9             11                 8                  9        ...
    Genome_2         2      1                  1           2             3                  2                  2        ...


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/hmm-hits-matrix-txt.md) to update this information.

