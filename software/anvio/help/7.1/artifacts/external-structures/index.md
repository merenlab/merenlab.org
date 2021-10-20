---
layout: page
title: external-structures [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/external-structures
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


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

By default, anvi'o predicts protein structures using MODELLER when creating a <span class="artifact-n">[structure-db](/software/anvio/help/7.1/artifacts/structure-db)</span>. Yet, if the user provides an external structures file, then anvi'o does not perform template-based homology modelling, and instead uses this file to obtain the structure information for the <span class="artifact-n">[structure-db](/software/anvio/help/7.1/artifacts/structure-db)</span>.

External structures is a user-provided TAB-delimited file that should follow this format:

|gene_callers_id|path|
|:---:|:---|
|1|path/to/gene1/structure.pdb|
|2|path/to/gene2/structure.pdb|
|3|path/to/gene3/structure.pdb|
|4|path/to/gene4/structure.pdb|
|7|path/to/gene5/structure.pdb|
|8|path/to/gene6/structure.pdb|
|(...)|(...)|

Each path should point to a <span class="artifact-n">[protein-structure-txt](/software/anvio/help/7.1/artifacts/protein-structure-txt)</span>.

{:.notice}
Please note that anvi'o will try its best to test the integrity of each file, and work with any limitations, however ultimately the user may be subject to the strict requirements set forth by anvi'o. For example, if a structure has a missing residue, you will hear about it.



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/external-structures.md) to update this information.

