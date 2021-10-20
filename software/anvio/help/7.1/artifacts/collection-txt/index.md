---
layout: page
title: collection-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/collection-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-collection](../../programs/anvi-export-collection)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-collection](../../programs/anvi-import-collection)</span> <span class="artifact-r">[anvi-script-get-coverage-from-bam](../../programs/anvi-script-get-coverage-from-bam)</span> <span class="artifact-r">[anvi-script-merge-collections](../../programs/anvi-script-merge-collections)</span></p>


## Description

This is a two-column TAB-delimited file without a header that describes a <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> by associating items with bin names. It can be used to import or export collections in and out of anvi'o databases, and/or transferring them between anvi'o projects. 

The first column in the file lists item names and the second column associates a given item with a bin. 

<div class="codeblock" markdown="1">
item_01    bin_1
item_02    bin_1
item_03    bin_1
item_04    bin_2
item_05    bin_3
item_06    bin_3
</div>

### The optinal bins info file

In addition to the essential file above, you can associate an optional TAB-delmited file with three columns with a collection to provide information about 'bins' in it, such as their source, and/or color to be used when they are displayed in <span class="artifact-n">[summary](/software/anvio/help/7.1/artifacts/summary)</span> outputs or anvi'o <span class="artifact-n">[interactive](/software/anvio/help/7.1/artifacts/interactive)</span> interfaces. Here is an example:

```
bin_1	N/A	 #c9d433
bin_2	N/A	 #e86548
bin_3	N/A	 #0b8500
```

In this file format the first column is a bin name, the second column is a source, and the third column is an HTML color.

{:.notice}
The source is a free form text and can be anything. We often use `anvi-interactive` or `CONCOCT` or `anvi-refine` for our bins to track which ones were manually refined, and which ones were coming from an automated binning algorithm.

You can provide this optional file to the program <span class="artifact-n">[anvi-import-collection](/software/anvio/help/7.1/programs/anvi-import-collection)</span> with the parameter `--bins-info`.

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/collection-txt.md) to update this information.

