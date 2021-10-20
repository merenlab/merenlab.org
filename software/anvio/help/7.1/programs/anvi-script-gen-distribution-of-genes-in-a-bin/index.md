---
layout: page
title: anvi-script-gen-distribution-of-genes-in-a-bin [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-gen-distribution-of-genes-in-a-bin
image:
  featurerelative: ../../../images/header.png
  display: true
---

Quantify the detection of genes in genomes in metagenomes to identify the environmental core. This is a helper script for anvi&#x27;o metapangenomic workflow.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[view-data](../../artifacts/view-data) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-items-txt](../../artifacts/misc-data-items-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program computes the detection of genes (inputted as a <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>) across your samples, so that you can visualize them in the <span class="artifact-n">[interactive](/software/anvio/help/7.1/artifacts/interactive)</span> interface. 

This program is used in [the metapangenomic workflow](https://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes) on genes with metagenomes as samples to visually identify the environmental core genes and accessory genes. 

### Inputs  

Essentially, you provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> pair, as well as the <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> you want to look at, and this program will  search each gene in your bin against the samples denoted in your <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>: 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen&#45;distribution&#45;of&#45;genes&#45;in&#45;a&#45;bin &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \ 
                                               &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                                               &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> \
                                               &#45;b <span class="artifact&#45;n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> 
</div>

There are two other parameters that you can set to focus the genes that you're looking at: 
- The minimum detection required for a gene to be included (by default, a gene must have a detection value of `0.5` in at least one of your samples)
-The minimum coverage required for a gene to be included (by default, a gene must have a total coverage of `0.25` times the mean total coverage in your data) 

### Outputs

This program will produce two outputs: 

1. `[your bin name]-GENE-COVs.txt`, which is a <span class="artifact-n">[view-data](/software/anvio/help/7.1/artifacts/view-data)</span> artifact. This is a matrix where each row represents a gene, each column represents one of your samples, and the cells each contain a coverage value. 
2. `[your bin name]-ENV-DETECTION.txt`, which is a <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span>. It is a two-column file, where each row is a gene and and the second column describes whether or not that gene is systematically detected in your samples. Thus, this can be added as an additional layer in the interface that describes describes which genes are detected in your samples. (as an example, see the outermost layer [here](https://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes))

Thus, after running this program on a bin with name `BIN_NAME`, you can run 

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;interactive](/software/anvio/help/7.1/programs/anvi&#45;interactive)</span> &#45;d BIN_NAME&#45;GENE&#45;COVs.txt \
                 &#45;A BIN_NAME&#45;ENV&#45;DETECTION.txt \
                 &#45;&#45;manual \
                 &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span>
</div>                                                   

This will visually show you the coverage and detection of your genes across your samples in the <span class="artifact-n">[interactive](/software/anvio/help/7.1/artifacts/interactive)</span> interface (simlarly to [this figure](https://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes)). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen-distribution-of-genes-in-a-bin.md) to update this information.


## Additional Resources


* [This program in action as part of the metapangenomic workflow](http://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-distribution-of-genes-in-a-bin) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
