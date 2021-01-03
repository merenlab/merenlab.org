---
layout: page
title: functions [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-import-functions](../../programs/anvi-import-functions)</span> <span class="artifact-p">[anvi-run-kegg-kofams](../../programs/anvi-run-kegg-kofams)</span> <span class="artifact-p">[anvi-run-ncbi-cogs](../../programs/anvi-run-ncbi-cogs)</span> <span class="artifact-p">[anvi-run-pfams](../../programs/anvi-run-pfams)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-analyze-synteny](../../programs/anvi-analyze-synteny)</span> <span class="artifact-r">[anvi-export-functions](../../programs/anvi-export-functions)</span></p>


## Description

This artifact **contains information about the functions of the genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>.**

It is a table within your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> that contains functional annotations for your genes. 

To get one of these for your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>, you can either import it (using <span class="artifact-n">[anvi-import-functions](/software/anvio/help/7/programs/anvi-import-functions)</span>) or make one yourself by running your contigs against one of two databases available in anvi'o:
* NCBI [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/) -- see <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/7/programs/anvi-run-ncbi-cogs)</span> for instructions
* EBI's [Pfam database](https://pfam.xfam.org/) -- see <span class="artifact-n">[anvi-run-pfams](/software/anvio/help/7/programs/anvi-run-pfams)</span> for instructions
* The [KOfamKOALA database](https://www.genome.jp/tools/kofamkoala/) -- see <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span> for instructions. Functions specifically from the KOfam database are used to run <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> and are also called <span class="artifact-n">[kegg-functions](/software/anvio/help/7/artifacts/kegg-functions)</span>.

You can also use EggNOG or InterProScan and then import the results into anvi'o, as described in [this blog post](http://merenlab.org/2016/06/18/importing-functions/).

This is used to run <span class="artifact-n">[anvi-analyze-synteny](/software/anvio/help/7/programs/anvi-analyze-synteny)</span>. 

You can also use <span class="artifact-n">[anvi-export-functions](/software/anvio/help/7/programs/anvi-export-functions)</span> to obtain a file containing these functional annotations through a <span class="artifact-n">[functions-txt](/software/anvio/help/7/artifacts/functions-txt)</span> artifact. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/functions.md) to update this information.

