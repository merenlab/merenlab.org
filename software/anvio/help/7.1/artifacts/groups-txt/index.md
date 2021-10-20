---
layout: page
title: groups-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/groups-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-compute-functional-enrichment-across-genomes](../../programs/anvi-compute-functional-enrichment-across-genomes)</span> <span class="artifact-r">[anvi-compute-metabolic-enrichment](../../programs/anvi-compute-metabolic-enrichment)</span> <span class="artifact-r">[anvi-display-functions](../../programs/anvi-display-functions)</span></p>


## Description

This is a 2-column TAB-delimited text file to associate a given set of items with a set of groups. Depending on the context, items here may be individual samples or genomes. The first column can have a header name `item`, `sample`, `genome` or anything else that is appropriate, and list the items that are relevant to your input data. The second column should have the header `group`, and associate each item in your data with a group.

Each item should be associated with a single group, and it is always a good idea to define groups using single words without any fancy characters. For instance, `HIGH_TEMPERATURE` or `LOW_FITNESS` are good group names. In contrast, `my group #1` or `IS-THIS-OK?`, are not quite good names for groups and may cause issues downstream depending on who uses this file.

Here is an example <span class="artifact-n">[groups-txt](/software/anvio/help/7.1/artifacts/groups-txt)</span> file:

|item|group|
|:--|:--|
|item_01|GROUP_A|
|item_02|GROUP_B|
|item_03|GROUP_A|
|(...)|(...)|

{:.warning}
If you are passing this file to the program <span class="artifact-n">[anvi-compute-metabolic-enrichment](/software/anvio/help/7.1/programs/anvi-compute-metabolic-enrichment)</span>, the names in the `sample` column must match those in the "modules" mode output file that you provide to the program via the `--modules-txt` parameter. If you know that the sample names match but you are still getting errors, you might need to specify which column in the "modules" mode output contains those sample names using the `--sample-header` parameter.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/groups-txt.md) to update this information.

