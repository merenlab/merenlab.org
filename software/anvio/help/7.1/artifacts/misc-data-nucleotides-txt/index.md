---
layout: page
title: misc-data-nucleotides-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/misc-data-nucleotides-txt
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

This a tab-delimited text file that describes information contained in a <span class="artifact-n">[misc-data-nucleotides](/software/anvio/help/7.1/artifacts/misc-data-nucleotides)</span>. 

To import this information into a database, use <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>. 

In this table, the first column should provide two pieces of information, both identifying a specific nucleotide position: the name of the contig the nucleotide is on, and its position on that contig. These should be separated by a colon. The following columns can contain any categorical or numerical data of your choosing.

Here is an example with very abstract data:

    item_name   categorical_data  numerical_data    data_group
    contig_1:4       group_1          4.3245         cool_data 
    contig_4:72      group_2          1.3542         cool_data
    contig_7:24      group_1          3.2526         cool_data
    ...

For a more concrete example, check out the example table for <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/7.1/artifacts/misc-data-amino-acids)</span> (which is formatted very similarly)  [here](http://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database). The second table on this page is what you would provide to <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-nucleotides-txt.md) to update this information.

