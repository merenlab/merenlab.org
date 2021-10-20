---
layout: page
title: internal-genomes [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/internal-genomes
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-gen-genomes-file](../../programs/anvi-script-gen-genomes-file)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-compute-functional-enrichment-across-genomes](../../programs/anvi-compute-functional-enrichment-across-genomes)</span> <span class="artifact-r">[anvi-compute-genome-similarity](../../programs/anvi-compute-genome-similarity)</span> <span class="artifact-r">[anvi-compute-metabolic-enrichment](../../programs/anvi-compute-metabolic-enrichment)</span> <span class="artifact-r">[anvi-dereplicate-genomes](../../programs/anvi-dereplicate-genomes)</span> <span class="artifact-r">[anvi-display-functions](../../programs/anvi-display-functions)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-gen-genomes-storage](../../programs/anvi-gen-genomes-storage)</span> <span class="artifact-r">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span> <span class="artifact-r">[anvi-meta-pan-genome](../../programs/anvi-meta-pan-genome)</span> <span class="artifact-r">[anvi-script-gen-functions-per-group-stats-output](../../programs/anvi-script-gen-functions-per-group-stats-output)</span> <span class="artifact-r">[anvi-script-gen-hmm-hits-matrix-across-genomes](../../programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span></p>


## Description

An internal genome is any <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> described in an anvi'o <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> stored in an anvi'o <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>. You can obtain one of these by binning a metagenome assembly (stored in an anvi'o <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>), which you can do either manually in the interactive interface or automatically with a binning software, and saving or importing it into a <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>.

The internal genomes file format enables anvi'o to work with one or more bins from one or more collections that may be defined in different anvi'o <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> files. A TAB-delimited internal genomes file will be composed of at least the following five columns:

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|Name_01|Bin_id_01|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_02|Bin_id_02|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_03|Bin_id_03|Collection_B|/path/to/another_profile.db|/path/to/another/contigs.db|
|(...)|(...)|(...)|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see **<span class="artifact-n">[external-genomes](/software/anvio/help/7.1/artifacts/external-genomes)</span>** and **<span class="artifact-n">[metagenomes](/software/anvio/help/7.1/artifacts/metagenomes)</span>**.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/internal-genomes.md) to update this information.

