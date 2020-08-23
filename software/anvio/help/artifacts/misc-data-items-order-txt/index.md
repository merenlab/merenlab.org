---
layout: post
title: misc-data-items-order-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-items-order](../../programs/anvi-export-items-order)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-items-order](../../programs/anvi-import-items-order)</span></p>

## Description

This is a text file that **contains the information in a <span class="artifact-n">[misc-data-items-order](/software/anvio/help/artifacts/misc-data-items-order)</span>**, used for importing into and exporting this information from your anvi'o project. 

If you open this file, it will look something like this: 

    (contig_4, ((contig_1, contig_2), contig_3))
    
    
(but probably much more complicated.) If this is imported into an Anvi'o project, the contigs will be displayed in the order `contig_4, contig_1, contig_2, contig_3`, and the following tree will be generated in the interface: 
    
    contig_4    contig_1    contig_2    contig_3
        |           |           |           |
        |           -------------           |
        |                 |                 |
        |                 -------------------
        |                           |
        -----------------------------
                    |
                    |


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-items-order-txt.md) to update this information.

