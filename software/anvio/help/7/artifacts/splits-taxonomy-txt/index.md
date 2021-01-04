---
layout: page
title: splits-taxonomy-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-splits-taxonomy](../../programs/anvi-export-splits-taxonomy)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This is the output of <span class="artifact-n">[anvi-export-splits-taxonomy](/software/anvio/help/7/programs/anvi-export-splits-taxonomy)</span>. It is in the same format as a <span class="artifact-n">[gene-taxonomy-txt](/software/anvio/help/7/artifacts/gene-taxonomy-txt)</span>, namely the first column identifies splits and the following columns describe the taxonomy hit associated with that split. 

For example:

    split_id  t_domain     t_phylum       t_class      ...
    1         Eukarya      Chordata       Mammalia
    2         Prokarya     Bacteroidetes  Bacteroidia
    ...



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/splits-taxonomy-txt.md) to update this information.

