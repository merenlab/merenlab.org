---
layout: page
title: misc-data-items-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span> <span class="artifact-p">[anvi-script-gen-distribution-of-genes-in-a-bin](../../programs/anvi-script-gen-distribution-of-genes-in-a-bin)</span> <span class="artifact-p">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-r">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Description

This a tab-delimited text file that describes information contained in a <span class="artifact-n">[misc-data-items](/software/anvio/help/7/artifacts/misc-data-items)</span>. 

To import this information into a database, use <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7/programs/anvi-import-misc-data)</span>. 

In this table, the first column should match the names of the items that you're displaying, and the following columns can contain any categorical or numerical data of your choosing.  (You can even be fancy and display data as a stacked bar graph.)

For an example, check out [the table on this page](http://merenlab.org/2017/12/11/additional-data-tables/#items-additional-data-table).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-items-txt.md) to update this information.

