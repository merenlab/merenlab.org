---
layout: page
title: internal-genomes [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-compute-functional-enrichment](../../programs/anvi-compute-functional-enrichment)</span> <span class="artifact-r">[anvi-compute-genome-similarity](../../programs/anvi-compute-genome-similarity)</span> <span class="artifact-r">[anvi-dereplicate-genomes](../../programs/anvi-dereplicate-genomes)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-gen-genomes-storage](../../programs/anvi-gen-genomes-storage)</span> <span class="artifact-r">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span> <span class="artifact-r">[anvi-meta-pan-genome](../../programs/anvi-meta-pan-genome)</span> <span class="artifact-r">[anvi-script-gen-hmm-hits-matrix-across-genomes](../../programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span></p>


## Description

An internal genome is any <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span> described in an anvi'o <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> stored in an anvi'o <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. You can obtain one of these by binning a metagenome assembly (stored in an anvi'o <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>), which you can do either manually in the interactive interface or automatically with a binning software, and saving or importing it into a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>.

The internal genomes file format enables anvi'o to work with one or more bins from one or more collections that may be defined in different anvi'o <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> files. A TAB-delimited internal genomes file will be composed of at least the following five columns:

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|Name_01|Bin_id_01|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_02|Bin_id_02|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_03|Bin_id_03|Collection_B|/path/to/another_profile.db|/path/to/another/contigs.db|
|(...)|(...)|(...)|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

## Additional columns

In some cases additional columns may be required to be in this file. Below is a table of the possible columns you may need.

| header | description | required for |
|----|----|----|
| group | name of the group that the genome belongs to (can be empty) | <span class="artifact-n">[anvi-compute-functional-enrichment](/software/anvio/help/7/programs/anvi-compute-functional-enrichment)</span> |

Also see **<span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span>** and **<span class="artifact-n">[metagenomes](/software/anvio/help/7/artifacts/metagenomes)</span>**.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/internal-genomes.md) to update this information.

