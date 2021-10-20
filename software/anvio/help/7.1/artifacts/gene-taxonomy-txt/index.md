---
layout: page
title: gene-taxonomy-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/gene-taxonomy-txt
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


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-taxonomy-for-genes](../../programs/anvi-import-taxonomy-for-genes)</span></p>


## Description

This is a file containing **the taxonomy information for the genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>**. 

You can use <span class="artifact-n">[anvi-import-taxonomy-for-genes](/software/anvio/help/7.1/programs/anvi-import-taxonomy-for-genes)</span> to integrate this information into your contigs database. See [this blog post](http://merenlab.org/2016/06/18/importing-taxonomy/) for a comprehensive tutorial. 

In its simplest form, this file is a tab-delimited text file that lists gene caller IDs and their associated taxonomy information. However, anvi'o can also parse outputs from taxonomy-based software like [Kaiju](https://github.com/bioinformatics-centre/kaiju) or [Centrifuge](https://github.com/infphilo/centrifuge). 

For example:

    gene_caller_id  t_domain     t_phylum       t_class      ...
          1         Eukarya      Chordata       Mammalia
          2         Prokarya     Bacteroidetes  Bacteroidia
          ...




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/gene-taxonomy-txt.md) to update this information.

