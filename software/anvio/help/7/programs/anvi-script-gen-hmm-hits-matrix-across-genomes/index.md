---
layout: page
title: anvi-script-gen-hmm-hits-matrix-across-genomes [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

A simple script to generate a TAB-delimited file that reports the frequency of HMM hits for a given HMM source across contigs databases.

See **[program help menu](../../../../vignette#anvi-script-gen-hmm-hits-matrix-across-genomes)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[hmm-hits-matrix-txt](../../artifacts/hmm-hits-matrix-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-source](../../artifacts/hmm-source) <img src="../../images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-hits](../../artifacts/hmm-hits) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program lets you look at the <span class="artifact-n">[hmm-hits](/software/anvio/help/7/artifacts/hmm-hits)</span> from a single <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span> across multiple genomes or bins, by creating a <span class="artifact-n">[hmm-hits-matrix-txt](/software/anvio/help/7/artifacts/hmm-hits-matrix-txt)</span>. 

The input of this program can be either an <span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span> or an <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span>. 

Here are two example run on an internal-genomes: 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen&#45;hmm&#45;hits&#45;matrix&#45;across&#45;genomes &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7/artifacts/internal&#45;genomes)</span> \
                                               &#45;&#45;hmm&#45;source Bacteria_71 \
                                               &#45;o output.txt
</div>

To list the <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span>s common to the datasets that you're analyzing, just add the flag `--list-hmm-sources`, as so: 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen&#45;hmm&#45;hits&#45;matrix&#45;across&#45;genomes &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/7/artifacts/external&#45;genomes)</span> \
                                               &#45;&#45;list&#45;hmm&#45;sources 
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen-hmm-hits-matrix-across-genomes.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-hmm-hits-matrix-across-genomes) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
