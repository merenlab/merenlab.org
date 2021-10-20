---
layout: page
title: anvi-plot-trnaseq [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-plot-trnaseq
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program to write plots of coverage and modification data from flexible groups of tRNA-seq seeds.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[trnaseq-contigs-db](../../artifacts/trnaseq-contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[trnaseq-seed-txt](../../artifacts/trnaseq-seed-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[modifications-txt](../../artifacts/modifications-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[trnaseq-plot](../../artifacts/trnaseq-plot) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **generates plots of groups of tRNA-seq seeds. The plots show seed coverages and the nucleotide frequencies at modification sites in each sample**.

The <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7.1/artifacts/trnaseq-workflow)</span> predicts tRNA seeds and their nucleotide modifications from a set of samples. The inspect webpage in <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> displays information on a selected seed, including coverages, mutation frequencies at predicted modification sites, and indel frequencies. anvi-plot-trnaseq generates plots that similarly show seed coverages in each sample but otherwise differ in many ways from the inspect page.

This program generates plots for a user-defined group of seeds. Seeds may be grouped by tRNA taxonomy and/or amino acid/anticodon identity. All of the seeds in a taxonomic/anticodon group are represented on the same plot. For example, all Arg-ACG seeds resolving to family Lachnospiraceae can be displayed on a single plot. Each panel of the plot shows coverages of all the seeds in a given sample. If there are five Lachnospiraceae Arg-ACG seeds, then five coverage traces will be stacked atop each other in each subplot, the seed with the highest mean coverage on the bottom. Nucleotide frequencies at predicted modification positions are shown as bars, with sections of the bar for each of the four nucleotides. If three of the five seeds have a modification, say m1A22, then the total height of the bar will rise to the height of the summed coverage of the three seeds at position 22. The number of seeds represented by the group is displayed on the plot, and the mean and 3' (discriminator nucleotide) coverages of the group of seeds is displayed adjacent to each sample subplot.

Multiple groups may be specified at the same time, producing a set of plots. If **only** a taxonomic group is given, then plots for **every** isoacceptor will be produced, e.g., Lachnospiraceae Arg-ACG, Arg-CCG, Arg-GCG, etc. If a taxonomic **rank** and anticodon is given, then plots for each taxon will be produced, e.g., at the family level, Lachnospiraceae Arg-ACG, Ruminococcaceae Arg-ACG, Bacteroidaceae Arg-ACG, etc.

anvi-plot-trnaseq is interactive through the command prompt, allowing the required <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span>, <span class="artifact-n">[seeds-specific-txt](/software/anvio/help/7.1/artifacts/seeds-specific-txt)</span>, and <span class="artifact-n">[modifications-txt](/software/anvio/help/7.1/artifacts/modifications-txt)</span> to be loaded only once and plots to be generated on the fly. Aesthetic parameters of the plots can be tweaked through the program. A comprehensive help menu with examples appears upon starting the program.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-plot-trnaseq.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-plot-trnaseq) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
