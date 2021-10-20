---
layout: page
title: misc-data-layer-orders-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/misc-data-layer-orders-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span></p>


## Description

This a tab-delimited text file that describes information contained in a <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7.1/artifacts/misc-data-layer-orders)</span>. 

To import this information into a database, use <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>. 

This table should contain trees formatted in either basic or newick form, where each branch represents the samples displayed by your layers. The order of the branches from left to right is the order they will be displayed in, from the center moving out. 

For an example, check out [the table on this page](http://merenlab.org/2017/12/11/additional-data-tables/#layer-orders-additional-data-table).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-layer-orders-txt.md) to update this information.

