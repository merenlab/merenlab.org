---
layout: page
title: view-data [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/view-data
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-gen-distribution-of-genes-in-a-bin](../../programs/anvi-script-gen-distribution-of-genes-in-a-bin)</span> <span class="artifact-p">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-matrix-to-newick](../../programs/anvi-matrix-to-newick)</span> <span class="artifact-r">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Description

View data refers to a matrx where each column represents a specific sample and each row describes some attribute of that sample (most often a sequence's abundance per sample). 

For example, in the [pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions), the `PROCHLORO-functions-occurrence-frequency.txt` is a view-data. 

You can use this to compute a distance matrix to generate a dendrogram (using <span class="artifact-n">[anvi-matrix-to-newick](/software/anvio/help/7.1/programs/anvi-matrix-to-newick)</span>) or direclty input it to <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> to visualize the distribution of your items across samples. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/view-data.md) to update this information.

