---
layout: page
title: anvi-gen-gene-consensus-sequences [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-gen-gene-consensus-sequences
image:
  featurerelative: ../../../images/header.png
  display: true
---

Collapse variability for a set of genes across samples.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[genes-fasta](../../artifacts/genes-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **provides consensus sequences for the genes within a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> pair**.

In other words, this collapses variability by assigning the most abundant nucleotide in your sample at each position, giving single consensus sequences for each gene for each sample. 

A basic run of this program will resemble the following: 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;gene&#45;consensus&#45;seuqences &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                                  &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                                  &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/7.1/artifacts/genes&#45;fasta)</span> 
</div>

The default output is a <span class="artifact-n">[genes-fasta](/software/anvio/help/7.1/artifacts/genes-fasta)</span>, but you can also get a tab-delimited output matrix by adding the flag  `--tab-delimited`.

You also have the option to focus on a subset of the data in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> by providing either: 

- A list of gene caller IDs (either as a parameter or through a file with one gene caller ID put line)
- A list of samples to focus on (as a file with a single sample name per line) 

### Additional Parameters 

- You have the option to change the variability engine (i.e. to codons), where variability at this level will be resolved. 
- To compress all variability profiles for each of your samples for a single gene, use the flag `--conpress samples`. This way, the program will only report one consensus sequence for each gene instead of reporting one for each sample. 
- You can get consensus sequences for each contig instead of for each gene with `--contigs-mode`
- To report all consensus sequences (even when there are no variable positions), activate `--quince-mode`


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-gene-consensus-sequences.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-gene-consensus-sequences) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
