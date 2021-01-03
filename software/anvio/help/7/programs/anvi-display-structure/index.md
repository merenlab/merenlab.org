---
layout: page
title: anvi-display-structure [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Interactively visualize sequence variants on protein structures.

See **[program help menu](../../../../vignette#anvi-display-structure)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[splits-txt](../../artifacts/splits-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage



This program opens an interactive interface to explore single amino acid variants (SAAVs) and single codon variants (SCVs) in the context of predicted tertiary protein structures and binding sites.  There are many example uses [here](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#display-metagenomic-sequence-variants-directly-on-predicted-structures) and you can work through an example as part of [the infant gut tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures) as well.  This is an integral program of anvi'o structure, which you can learn more about [here](https://merenlab.org/software/anvio-structure/).


In short, this program enables users to explore sequence variation in the context of 3D protein structure, which reveals insight that cannot be learned from purely sequence-based approaches.


### Before running 

To run this program, you'll need to have created a <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span> which can be easily done with a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> and the program <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/7/programs/anvi-gen-structure-database)</span>.


You'll also need a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> that was created using <span class="artifact-n">[anvi-profile](/software/anvio/help/7/programs/anvi-profile)</span>'s flag `--profile-SCVs`, which means that single codon variants (SCVs) have been profiled. Very sorry if this forces you to re-profile, but as of v6.2, this is now a very expedient process.


### Basic Run 

There are two ways to provide the variability information to this program.  

The first is to provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> pair, and let this program calculate SAAVs and SCVs as they are requested by the interface.


<div class="codeblock" markdown="1">
anvi&#45;display&#45;structure &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/7/artifacts/structure&#45;db)</span> \
                       &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                       &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> 
</div>

The second is to use <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/7/programs/anvi-gen-variability-profile)</span> to create a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7/artifacts/variability-profile-txt)</span>. This way, you pre-load all of the variability data and don't have to wait for <span class="artifact-n">[anvi-display-structure](/software/anvio/help/7/programs/anvi-display-structure)</span> to calculate variability on-the-fly. This option is probably most convenient in instances where you have already generated a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7/artifacts/variability-profile-txt)</span> for other reasons. If you fall into this camp, you can run <span class="artifact-n">[anvi-display-structure](/software/anvio/help/7/programs/anvi-display-structure)</span> as so:


<div class="codeblock" markdown="1">
anvi&#45;display&#45;structure &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/7/artifacts/structure&#45;db)</span> \
                       &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                       &#45;v <span class="artifact&#45;n">[variability&#45;profile&#45;txt](/software/anvio/help/7/artifacts/variability&#45;profile&#45;txt)</span>
</div>

{:.notice}
You still must provide the <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> used to generate the <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7/artifacts/variability-profile-txt)</span>, since it contains other necessary information such as functional annotations and ligand binding predictions.  You may optionally provide a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> if custom sample grouping is important to you.

{:.notice}
During <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/7/programs/anvi-gen-variability-profile)</span>, if you are _only_ interested in genes that have predicted structures, you may want to run <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/7/programs/anvi-gen-variability-profile)</span> with the flag `--only-if-structure`.

### Refining your search

You have several options to refine what proteins and variants you're looking at: 

- Provide a list of gene caller IDs to only display specific genes (this can be provided either directly as a parameter or as a file with one gene caller ID per line)
- Specify the minimum departure from the consensus sequence. This is a number from 0-1 that describes the threshold for a variability position to be displayed. For example, if this is set to 0.2, then all SAAVs and SCVs where less than 20 percent of the reads vary from the consensus sequence will not be displayed.
- Specify samples of interest. Those in your <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> or <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7/artifacts/variability-profile-txt)</span> that are not in the samples of interest will be filtered out.

If you're choosing to have <span class="artifact-n">[anvi-display-structure](/software/anvio/help/7/programs/anvi-display-structure)</span> calculate variability on-the-fly, you can speed things up by choosing to _only_ calculate SAAVs or _only_ calculate SCVs.


### Other parameters 

Power users can also change the server configuration (i.e. set the IP address, port number, browser path, server password, etc.)




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-display-structure.md) to update this information.


## Additional Resources


* [The overview page from the release](http://merenlab.org/software/anvio-structure/)

* [The section of the Infant Gut Tutorial focused on anvi-display-structure](http://merenlab.org/tutorials/infant-gut/#chapter-vii-from-single-amino-acid-variants-to-protein-structures)

* [Integrating sequence variants and predicted protein structures](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-display-structure) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
