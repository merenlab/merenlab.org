---
layout: page
title: anvi-gen-fixation-index-matrix [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-gen-fixation-index-matrix
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate a pairwise matrix of a fixation indices between samples.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[splits-txt](../../artifacts/splits-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[fixation-index-matrix](../../artifacts/fixation-index-matrix) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage



This program generates a matrix of the pairwise fixation indices (F<sub>ST</sub>) between your samples.

### What's a fixation index?

As described [in the Infant Gut Tutorial](https://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst), the fixation index is a measure of the distance between two populations, based on their sequence variants (usually SNVs). Specifically, the fixation index is the ratio between the variance in allele frequency between subpopulations and the variance in the total population. 


The fixation index has its own [Wikipedia page](https://en.wikipedia.org/wiki/Fixation_index) and is a special case of [F-statistics](https://en.wikipedia.org/wiki/F-statistics). 


In anvi'o, the fixation index is calculated in accordance with [Schloissnig et al.  (2013)](https://doi.org/10.1038/nature11711)'s work to allow variant positions with multiple competing alleles.


## Anvi-gen-fixation-index 

There are two ways to run this program.  

### Input 1: Variability Profile 

The simplest one is the one shown [in the Infant Gut Tutorial](https://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst): just provide a <span class="artifact-n">[variability-profile](/software/anvio/help/7.1/artifacts/variability-profile)</span>, like so: 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;fixation&#45;index&#45;matrix &#45;&#45;variability&#45;profile <span class="artifact&#45;n">[variability&#45;profile](/software/anvio/help/7.1/artifacts/variability&#45;profile)</span> \
                               &#45;&#45;output&#45;file my_matrix.txt
</div>

This will use the information in your <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7.1/artifacts/variability-profile-txt)</span> to generate the fixation index for each of the pairwise sample comparisons, and store the results in a <span class="artifact-n">[fixation-index-matrix](/software/anvio/help/7.1/artifacts/fixation-index-matrix)</span> named `my_matrix.txt`.  

### Input 2: Anvi'o databases

Instead of providing a <span class="artifact-n">[variability-profile](/software/anvio/help/7.1/artifacts/variability-profile)</span>, you can instead provide the inputs to <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/7.1/programs/anvi-gen-variability-profile)</span> and let anvi'o do all of the work for you. Specifically, this means providing a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> pair to find your variability positions and a specific subset to focus on in any of these ways: 

- Provide a list of gene caller IDs (as a parameter with the flag `--gene-caller-ids` or in a file with one ID per line with the flag `--genes-of-interest`)
- Provide a list of splits (in a <span class="artifact-n">[splits-txt](/software/anvio/help/7.1/artifacts/splits-txt)</span>)
- Provide a <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> and <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>

Additionally, you can add structural annotations by inputting a <span class="artifact-n">[structure-db](/software/anvio/help/7.1/artifacts/structure-db)</span> (and focus only on genes with structural annotations with the flag `--only-if-structure`) or choose to focus on only a subset of your samples by providing a file of samples of interest.  

When doing this, you can also set the variability engine to get the fixation index for SCVs (`--engine CDN`) or SAAVs (`--engine AA`). 

You can find more information about these parameters on the page for <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/7.1/programs/anvi-gen-variability-profile)</span>. 

### Additional Parameters

While a fixation index is usually between 0 and 1, it is possible for an index to be negative (usually because of out-breeding). By default, anvi'o sets these negative values to 0, but you can choose to keep the negative values with the flag `--keep-negatives` 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-fixation-index-matrix.md) to update this information.


## Additional Resources


* [Utilizing fixation index to study SAR11 population structure](http://merenlab.org/data/sar11-saavs/#generating-distance-matrices-from-fixation-index-for-saavs-and-snvs-data)

* [Measuring Distances Between Genomes in the Infant Gut Tutorial](http://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-fixation-index-matrix) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
