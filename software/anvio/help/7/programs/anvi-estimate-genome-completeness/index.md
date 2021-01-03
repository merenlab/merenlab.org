---
layout: page
title: anvi-estimate-genome-completeness [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Estimate completion and redundancy using domain-specific single-copy core genes.

See **[program help menu](../../../../vignette#anvi-estimate-genome-completeness)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[completion](../../artifacts/completion) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>

## Usage


This program estimates the completeness and redundancy of genomes provided to it, based on domain-level single-copy core genes. 

{:.notice}
Wondering what single-copy core genes anvi'o uses? Check out <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span>. It uses the tables populated when you ran <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7/programs/anvi-run-hmms)</span> on your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>. 

Genomes provided to this program must be contained in either a <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span> (within a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>) or a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> (which can be provided alone or as part of an <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span>). 

### Running on contigs databases 

For example, calling 

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;genome&#45;completeness &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> 
</div>

will output to the terminal the completition and redundancy of the single-copy core genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>, assuming that all of its contigs represent a single genome. To output this information to a file, you can add the flag `-o` and provide an output path. 

To get this information for several contigs databases at once, you can provide them as an <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span>, as so:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;genome&#45;completeness &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/7/artifacts/external&#45;genomes)</span> \
                                  &#45;o completition.txt
</div>

### Running on bins 

To get this data for a series of bins, just provide a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> and <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>. 

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;genome&#45;completeness &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                                  &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                                  &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7/artifacts/collection)</span> 
</div>

To see what collections are contained in your contigs database, call 

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;genome&#45;completeness &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                                  &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                                  &#45;&#45;list&#45;collections
</div>

or run <span class="artifact-n">[anvi-show-collections-and-bins](/software/anvio/help/7/programs/anvi-show-collections-and-bins)</span> for a more comprehensive overview. 

If you're looking for a more comprehensive overview of your entire collection and its contents, the completition and redunduncy statistics for your bins are also included when you run <span class="artifact-n">[anvi-summarize](/software/anvio/help/7/programs/anvi-summarize)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-estimate-genome-completeness.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-estimate-genome-completeness) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
