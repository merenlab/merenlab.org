---
layout: page
title: anvi-gen-genomes-storage [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-gen-genomes-storage
image:
  featurerelative: ../../../images/header.png
  display: true
---

Create a genome storage from internal and/or external genomes for a pangenome analysis.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **generates a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span>, which stores information about your genomes, primarily for use in pangenomic analysis.** 

Genomes storage databases are to Anvi'o's pangenomic workflow what a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> is to a metagenomic workflow: it stores vital information and is passed to most programs you'll want to run. 

Once you've generated a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span>, you can run <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/7.1/programs/anvi-pan-genome)</span>, which creates a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> and runs various pangenomic analyses (including calculating the similarities between your sequences, identifying gene clusters, and organizing your gene clusters and genomes). After that, you can display your pangenome with <span class="artifact-n">[anvi-display-pan](/software/anvio/help/7.1/programs/anvi-display-pan)</span> For more information, check out [the pangenomic workflow](http://merenlab.org/2016/11/08/pangenomics-v2/#generating-an-anvio-genomes-storage).

### Inputs: internal and external genomes

You can initialize your genomes storage database with <span class="artifact-n">[internal-genomes](/software/anvio/help/7.1/artifacts/internal-genomes)</span>, <span class="artifact-n">[external-genomes](/software/anvio/help/7.1/artifacts/external-genomes)</span>, or both. 

<span class="artifact-n">[internal-genomes](/software/anvio/help/7.1/artifacts/internal-genomes)</span> describe genomes that are described by a <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> within a <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> that is already within an Anvi'o <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>. For example, if you had gone through [the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/) and had several MAGs that you wanted to run pangenomic analyses on. 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;genomes&#45;storage &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7.1/artifacts/internal&#45;genomes)</span> \
                         &#45;o <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span>
</div>

{:.notice}
The name of your genomes storage database (which follows the `-o` flag) must end with `-GENOMES.db`. This just helps differenciate it from other types of Anvi'o databases, such as the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>. 

In contrast, <span class="artifact-n">[external-genomes](/software/anvio/help/7.1/artifacts/external-genomes)</span> describe genomes that are contained in a <span class="artifact-n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> file that you've turned into a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> (using <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span>).  For example, if you had downloaded genomes from [NCBI](https://www.ncbi.nlm.nih.gov/). 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;genomes&#45;storage &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/7.1/artifacts/external&#45;genomes)</span> \
                         &#45;o <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span>
</div>

You can also create a genomes storage database from both types of genomes at the same time. For example, if you had MAGs from a metagenomic analysis on an environmental sample and wanted to compare them with the reference genomes on [NCBI](https://www.ncbi.nlm.nih.gov/). To run this, simply provide both types of genomes as parameters, as so: 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;genomes&#45;storage &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7.1/artifacts/internal&#45;genomes)</span> \
                         &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/7.1/artifacts/external&#45;genomes)</span> \
                         &#45;o <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span>
</div>

### Changing the gene caller

By default, Anvi'o will use [Prodigal](https://github.com/hyattpd/Prodigal) and will let you know if you have gene calls identified by other gene callers. However, you are welcome to explicitly use a specific gene caller with the flag `--gene-caller`. 

If you're wondering what gene callers are available in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, you can check by running the program <span class="artifact-n">[anvi-export-gene-calls](/software/anvio/help/7.1/programs/anvi-export-gene-calls)</span> on a specific <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> with the flag `--list-gene-callers`. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-genomes-storage.md) to update this information.


## Additional Resources


* [A tutorial on pangenomics](http://merenlab.org/2016/11/08/pangenomics-v2/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-genomes-storage) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
