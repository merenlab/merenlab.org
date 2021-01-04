---
layout: page
title: misc-data-amino-acids [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-p">[anvi-run-interacdome](../../programs/anvi-run-interacdome)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-delete-misc-data](../../programs/anvi-delete-misc-data)</span> <span class="artifact-r">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span></p>


## Description


This is a section of your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> that contains custom additional information about specific amino acid residues.  

Take a look at [this blogpost](http://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database) for potential uses in the InteracDome (which will likely be added to anvi'o in v7) and the motivation behind this program.  

Similarly to other types of miscellaneous data (like <span class="artifact-n">[misc-data-items](/software/anvio/help/7/artifacts/misc-data-items)</span>), this information is either numerical or categorical and can be populated into a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> (from a <span class="artifact-n">[misc-data-amino-acids-txt](/software/anvio/help/7/artifacts/misc-data-amino-acids-txt)</span>) with <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7/programs/anvi-import-misc-data)</span>. It is also displayed when you run <span class="artifact-n">[anvi-show-misc-data](/software/anvio/help/7/programs/anvi-show-misc-data)</span> and can be exported or deleted with <span class="artifact-n">[anvi-export-misc-data](/software/anvio/help/7/programs/anvi-export-misc-data)</span> and <span class="artifact-n">[anvi-delete-misc-data](/software/anvio/help/7/programs/anvi-delete-misc-data)</span> respectively.  

For example, this could describe various key residues for binding to ligands, or residues otherwise determined to be important to the user for whatever reason.  



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-amino-acids.md) to update this information.

