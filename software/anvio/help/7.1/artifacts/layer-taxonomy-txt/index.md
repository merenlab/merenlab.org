---
layout: page
title: layer-taxonomy-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/layer-taxonomy-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-taxonomy-for-layers](../../programs/anvi-import-taxonomy-for-layers)</span></p>


## Description

This is a text file containing taxonomy information for your layers (the same information as a <span class="artifact-n">[layer-taxonomy](/software/anvio/help/7.1/artifacts/layer-taxonomy)</span>). You can bring this information into a <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span> using <span class="artifact-n">[anvi-import-taxonomy-for-layers](/software/anvio/help/7.1/programs/anvi-import-taxonomy-for-layers)</span>. 

This is a tab-delimited text file that is formatted similarly to a <span class="artifact-n">[gene-taxonomy-txt](/software/anvio/help/7.1/artifacts/gene-taxonomy-txt)</span>. The first column describes the names of your layers, and the following columns each correspond to the taxonomy level described in the header. Here is an example:

    sample  t_domain    t_phylum    t_class     ...
     c1     Eukaryea    Chordata    Mammalia
     ...
     



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/layer-taxonomy-txt.md) to update this information.

