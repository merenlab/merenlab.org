---
layout: page
title: genbank-file [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/genbank-file
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-script-process-genbank](../../programs/anvi-script-process-genbank)</span></p>


## Description

The GenBank file format was created by NCBI. 

You can find an [explination](https://www.ncbi.nlm.nih.gov/genbank/) and [example](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) on the NCBI website. 

In anvi'o, this is used by <span class="artifact-n">[anvi-script-process-genbank](/software/anvio/help/7.1/programs/anvi-script-process-genbank)</span> to convert the information in the genbank file to a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span>, <span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span>, and <span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/genbank-file.md) to update this information.

