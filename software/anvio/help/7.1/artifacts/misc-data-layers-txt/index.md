---
layout: page
title: misc-data-layers-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /m/misc-data-layers-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span> <span class="artifact-p">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-r">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Description

This a tab-delimited text file that describes information contained in a <span class="artifact-n">[misc-data-layers](/software/anvio/help/main/artifacts/misc-data-layers)</span>. 

To import this information into a database, use <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/main/programs/anvi-import-misc-data)</span>. 

In this table, the first column should match the names of the samples that you're displaying, and the following columns can contain any categorical or numerical data of your choosing. (You can even be fancy and display data as a stacked bar graph.)

For an example, check out [the table on this page](http://merenlab.org/2017/12/11/additional-data-tables/#layers-additional-data-table).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-layers-txt.md) to update this information.

