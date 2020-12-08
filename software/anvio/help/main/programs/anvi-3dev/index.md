---
layout: page
title: anvi-3dev [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Interactively visualize sequence variants on protein structures.

See **[program help menu](../../../vignette#anvi-3dev)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[splits-txt](../../artifacts/splits-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program opens the interactive interface to visualize variable positions directly on the 3D structure of a protein. There are many example uses [here](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#display-metagenomic-sequence-variants-directly-on-predicted-structures) and you can work through an example as part of [the infant gut tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures). 

In short, though, this interface lets you view the predicted 3D structures of your protein (as stored in a <span class="artifact-n">[structure-db](/software/anvio/help/main/artifacts/structure-db)</span>) with the variability positions (SCVs and SAAVs, usually determined from your metagenomic data) directly mapped on. This can give you new insights, for example the solvent accessibility of individual SAAVs and the strucutral distribution of variability positions. This is espeically useful after using <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/main/programs/anvi-import-misc-data)</span> to annotate additional data into your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, such as binding sites for substrates or other enzymes. 

### Before running

To run this program, you'll need to have used your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> created a <span class="artifact-n">[structure-db](/software/anvio/help/main/artifacts/structure-db)</span> with <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/main/programs/anvi-gen-structure-database)</span>. 

You'll also need a <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> that has had its SCVs profiled. In other words, when you ran either <span class="artifact-n">[anvi-profile](/software/anvio/help/main/programs/anvi-profile)</span> or <span class="artifact-n">[anvi-merge](/software/anvio/help/main/programs/anvi-merge)</span>, you need to have added the flag `--profile-SCVs` for this to work. This may be computationally intensive, but it is also necessary to run anvi-3dev. 

### Basic Run

There are two ways to provide the variability information to this program. 

The first is to provide a <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> pair, and let this program calculate the variability positions for you in the moment. 

<div class="codeblock" markdown="1">
anvi&#45;3dev &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/main/artifacts/structure&#45;db)</span> \
          &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
          &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> 
</div>

The second is to use <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/main/programs/anvi-gen-variability-profile)</span> to create a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/main/artifacts/variability-profile-txt)</span>. This way, you don't have to wait for <span class="artifact-n">[anvi-3dev](/software/anvio/help/main/programs/anvi-3dev)</span> to run all of this analysis and it is easier to share your variaiblity information. 

For this purpose, you'll probably want to run <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/main/programs/anvi-gen-variability-profile)</span> with the flag `--only-if-structure` so that it only calculates varaibility proteins that can be visualized. Then you can run <span class="artifact-n">[anvi-3dev](/software/anvio/help/main/programs/anvi-3dev)</span> as so:

<div class="codeblock" markdown="1">
anvi&#45;3dev &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/main/artifacts/structure&#45;db)</span> \
          &#45;v <span class="artifact&#45;n">[variability&#45;profile&#45;txt](/software/anvio/help/main/artifacts/variability&#45;profile&#45;txt)</span>
</div>

### Refining your search 

You have several options to refine what proteins and variable positions you're looking at: 

- Provide a list of gene caller IDs to only display specific genes (this can be provided either directly as a parameter or as a file with one gene caller ID per line)
- Provide a <span class="artifact-n">[splits-txt](/software/anvio/help/main/artifacts/splits-txt)</span> to only look at specfic splits (though usually you'll want to refine your search by genes)
- Specify the minimum departure from the consensus sequence. This is a number from 0-1 that describes the threshold for a variability position to be displayed. For example, if this is set to 0.2, then all SAAVs and SCVs where less than 20 percent of the reads vary from the consensus sequence will not be displayed. 

If you're choosing to not use a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/main/artifacts/variability-profile-txt)</span> and have anvi-3dev calculate your variability, you can also change a few other parameters. If you are using a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/main/artifacts/variability-profile-txt)</span>, you can instead set these parameters when running <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/main/programs/anvi-gen-variability-profile)</span>
- Provide a list of samples to use to calculate variable positions 
- Specify whether to look speicfically at SCVs or SAAVs. 

### Other parameters

Power users can also change the server configuration (i.e. set the IP address, port number, browser path, server password, etc.)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-3dev.md) to update this information.


## Additional Resources


* [The overview page from the release](http://merenlab.org/software/anvi-3dev/)

* [The section of the Infant Gut Tutorial focused on anvi-3dev](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures)

* [Integrating sequence variants and predicted protein structures](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-3dev) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
