---
layout: page
title: anvi-tabulate-trnaseq [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-tabulate-trnaseq
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program to write standardized tab-delimited files of tRNA-seq seed coverage and modification results.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[trnaseq-contigs-db](../../artifacts/trnaseq-contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[trnaseq-profile-db](../../artifacts/trnaseq-profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[trnaseq-seed-txt](../../artifacts/trnaseq-seed-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[modifications-txt](../../artifacts/modifications-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **generates tabular files of tRNA-seq seed coverage and modification data that are easily manipulable by the user**.

anvi-tabulate-trnaseq is part of the <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7.1/artifacts/trnaseq-workflow)</span>, and is run following the finalization of tRNA seeds by <span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span>.

This program generates a table, <span class="artifact-n">[seeds-specific-txt](/software/anvio/help/7.1/artifacts/seeds-specific-txt)</span>, containing the specific coverage of each nucleotide position in each seed in every sample. If a nonspecific <span class="artifact-n">[trnaseq-profile-db](/software/anvio/help/7.1/artifacts/trnaseq-profile-db)</span> is also provided, this program generates a table of nonspecific coverages, <span class="artifact-n">[seeds-non-specific-txt](/software/anvio/help/7.1/artifacts/seeds-non-specific-txt)</span>. The distinction between specific and nonspecific coverage is explained in the <span class="artifact-n">[trnaseq-profile-db](/software/anvio/help/7.1/artifacts/trnaseq-profile-db)</span> artifact. These coverage tables have one row per seed per sample. They have three header rows for different ways of describing tRNA nucleotide positions: canonical position name (e.g., "discriminator_1"), canonical position (e.g., "73"), and "ordinal" position relative to all the other **possible** positions (e.g., "95").

anvi-tabulate-trnaseq also generates a table, <span class="artifact-n">[modifications-txt](/software/anvio/help/7.1/artifacts/modifications-txt)</span>, containing information on each predicted modification position in each seed, with one row per modification per seed per sample. This table includes four columns of position coverage counts of the four nucleotides.

All tables include taxonomic annotations of the seeds; annotations are added to the <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span> by <span class="artifact-n">[anvi-run-trna-taxonomy](/software/anvio/help/7.1/programs/anvi-run-trna-taxonomy)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-tabulate-trnaseq.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-tabulate-trnaseq) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
