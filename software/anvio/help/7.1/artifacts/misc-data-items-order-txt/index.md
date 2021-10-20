---
layout: page
title: misc-data-items-order-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/misc-data-items-order-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-items-order](../../programs/anvi-export-items-order)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-items-order](../../programs/anvi-import-items-order)</span></p>


## Description

This is a text file that **contains the information for a <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7.1/artifacts/misc-data-items-order)</span>**, used for importing into and exporting this information from your anvi'o project.

## NEWICK order

If you intend to import a tree order, the contents of your file should look something like this (but probably much more complicated depending on the number of items in your anvi'o database): 

```
(contig_4, ((contig_1, contig_2), contig_3))
```

When a NEWICK order is imported into an anvi'o project, the contigs will be displayed in the order `contig_4, contig_1, contig_2, contig_3`, and the following tree will be generated in the interface:

```
    contig_4    contig_1    contig_2    contig_3
        |           |           |           |
        |           -------------           |
        |                 |                 |
        |                 -------------------
        |                           |
        -----------------------------
                    |
                    |
```

## LIST order

Alternative to the NEWICK order, you can provide a list of items in flat form. For instance, if you want to order your items this way, your text file should look like the following, where each line contains a single item name in your database:

```
contig_4
contig_1
contig_2
contig_3
```

{:.warning}
After importing an order into a database, you may need to specifically select that order in the interactive interface through the "Item orders" dropbox and re-draw your display to change the default order.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/misc-data-items-order-txt.md) to update this information.

