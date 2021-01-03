---
layout: page
title: anvi-update-structure-database [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Add or re-run genes from an already existing structure database. All settings used to generate your database will be used in this program.

See **[program help menu](../../../../vignette#anvi-update-structure-database)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program is used to add additional genes to or re-run the analysis of genes already within a <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span>.

For that reason, it is very similar to <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/7/programs/anvi-gen-structure-database)</span> and the parameters used to run that program (when you first generated your <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span>) will be automatically applied when you run this program. To know what MODELLER parameters are being used, you run this program on a <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span> with the flag `--list-modeller-params`. 

To run this program, just provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> and <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span>, and name your genes of interest (either in a file or directly). If the named genes are not already in your <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span>, they will be added to the database. 

For example, if your <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span> already contains the genes with caller-IDs 1, 2 and 3, and you run

<div class="codeblock" markdown="1">
anvi&#45;update&#45;structure&#45;database &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                               &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/7/artifacts/structure&#45;db)</span> \
                               &#45;&#45;gene&#45;caller&#45;ids 1,4,5
</div>

Then the structural analysis for genes 4 and 5 will be added to your <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span> (assuming templates are found). Gene 1 will be ignored, since it is already present.

If instead you want to re-run the structural analysis on genes that are already in your <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span>, you'll need to specify that by adding the flag `--rerun-genes`

<div class="codeblock" markdown="1">
anvi&#45;update&#45;structure&#45;database &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                               &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/7/artifacts/structure&#45;db)</span> \
                               &#45;&#45;gene&#45;caller&#45;ids 1,4,5 \
                               &#45;&#45;rerun&#45;genes
</div>

Now, the program will rerun the analysis for gene 1 and will still add genes 4 and 5 to the <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span>. 

Both of these runs will have the same MODELLER parameters as your run of <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/7/programs/anvi-gen-structure-database)</span>. However, to get the raw outputs, you will need to use the parameter `--dump-dir`. You can also set a specific MODELLER program with `--modeller-executable`. Parameters for multi-threading would also have to be given again.



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-update-structure-database.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-update-structure-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
