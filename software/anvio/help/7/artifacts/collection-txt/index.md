---
layout: page
title: collection-txt [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-collection](../../programs/anvi-export-collection)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-collection](../../programs/anvi-import-collection)</span> <span class="artifact-r">[anvi-script-get-coverage-from-bam](../../programs/anvi-script-get-coverage-from-bam)</span> <span class="artifact-r">[anvi-script-merge-collections](../../programs/anvi-script-merge-collections)</span></p>


## Description

This file describes the contents of a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> in a readable format. This is used for taking collections in and out of Anvi'o, or transferring them between Anvi'o projects. 

This is a tab-delimited file that contains two columns. The first lists all of the contigs contained within the bins in this collection. The second lists the associated bin that that contig is placed in. For examples, check out [this page](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_contigs.txt) or below.

<div class="codeblock" markdown="1">
contig_1    bin_1
contig_2    bin_1
contig_3    bin_1
contig_4    bin_2
contig_5    bin_3
contig_6    bin_3
</div>

You can also list splits instead of contigs in the left column, as seen [here](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_splits.txt).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/collection-txt.md) to update this information.

