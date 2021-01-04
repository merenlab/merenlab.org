---
layout: page
title: splits-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-compute-completeness](../../programs/anvi-compute-completeness)</span> <span class="artifact-r">[anvi-display-structure](../../programs/anvi-display-structure)</span> <span class="artifact-r">[anvi-gen-fixation-index-matrix](../../programs/anvi-gen-fixation-index-matrix)</span> <span class="artifact-r">[anvi-gen-variability-profile](../../programs/anvi-gen-variability-profile)</span> <span class="artifact-r">[anvi-get-aa-counts](../../programs/anvi-get-aa-counts)</span></p>


## Description

This is a **text file containing a list of splits**. 

This file has only one column with one split name per line. For example 

    split_name_1
    split_name_2
    split_name_3
    
    
This kind of file is used when you want to focus your analysis on only a specific set of splits. Just provide one of these and the program will only look at `split_name_1`, `split_name_2`, and `split_name_3`.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/splits-txt.md) to update this information.

