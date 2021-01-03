---
layout: page
title: anvi-compute-gene-cluster-homogeneity [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Compute homogeneity for gene clusters.

See **[program help menu](../../../../vignette#anvi-compute-gene-cluster-homogeneity)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program **computes both the geometric homogeneity and functional homogeneity for the gene clusters in a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span>.** 

*Geometric homogeneity* and *functional homogeneity* are anvi'o specific terms that describe how similar genes within a gene cluster are to each other in different ways. Briefly, geometric homogeneity compares the positions of gaps in the aligned residues without considering specific amino acids, and functional homogeneity examines point mutations to amino acids and compares how similar the resulting amino acids are chemically. See [this page](http://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters) for more details. 

You can run this program as so: 

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;gene&#45;cluster&#45;homogeneity &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span> \
                                      &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7/artifacts/genomes&#45;storage&#45;db)</span> \
                                      &#45;o path/to/output.txt \
                                      &#45;&#45;store&#45;in&#45;db
</div>

This run will put the output directly in the database, as well as provide it as a separate file as the specified output path. 

You also have the option to calculate this information about only specific gene clusters, either by providing a gene cluster ID, list of gene cluster IDs, <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> or <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span>. 

To save on runtime, you can also enable `--quick-homogeneity`, which will not check for horizontal geometric homogenity (i.e. it will not look at alignments within a single gene). This will be less accurate for detailed analyses, but it will run faster. 

Here is an example run that uses this flag and only looks at a specific collection: 

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;gene&#45;cluster&#45;homogeneity &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span> \
                                      &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7/artifacts/genomes&#45;storage&#45;db)</span> \
                                      &#45;o path/to/output.txt \
                                      &#45;&#45;store&#45;in&#45;db \ 
                                      &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7/artifacts/collection)</span> \
                                      &#45;&#45;quick&#45;homogeneity 
</div>

You can also use multithreading if you're familiar with that. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-compute-gene-cluster-homogeneity.md) to update this information.


## Additional Resources


* [The role of gene cluster homogeneity described in the Anvi&#x27;o pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-compute-gene-cluster-homogeneity) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
