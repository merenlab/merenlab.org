---
layout: page
title: anvi-get-sequences-for-gene-clusters [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-get-sequences-for-gene-clusters
image:
  featurerelative: ../../../images/header.png
  display: true
---

Do cool stuff with gene clusters in anvi&#x27;o pan genomes.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[genes-fasta](../../artifacts/genes-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[concatenated-gene-alignment-fasta](../../artifacts/concatenated-gene-alignment-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-items](../../artifacts/misc-data-items) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This aptly-named program **gets the sequences for the gene clusters stored in a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> and returns them as either a <span class="artifact-n">[genes-fasta](/software/anvio/help/7.1/artifacts/genes-fasta)</span> or a <span class="artifact-n">[concatenated-gene-alignment-fasta](/software/anvio/help/7.1/artifacts/concatenated-gene-alignment-fasta)</span>** (which you can use to run <span class="artifact-n">[anvi-gen-phylogenomic-tree](/software/anvio/help/7.1/programs/anvi-gen-phylogenomic-tree)</span>). This gives you advanced access to your gene clusters, which you can take out of anvi'o, use for phylogenomic analyses, or do whatever you please with. 

You also have the option to output the sequences of your choice as a <span class="artifact-n">[misc-data-items](/software/anvio/help/7.1/artifacts/misc-data-items)</span> (with `add-into-items-additional-data-table`), which can be added to the <span class="artifact-n">[interactive](/software/anvio/help/7.1/artifacts/interactive)</span> interface as additional layers. 

While the number of parameters may seem daunting, many of the options just help you specify exactly which gene clusters you want to get the sequences  from. 

### Running on all gene clusters

Here is a basic run, that will  export alignments for every single gene cluster found in the <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> as amino acid sequences :

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;clusters &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                     &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
                                     &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/7.1/artifacts/genes&#45;fasta)</span>
</div>

To get the DNA sequences instead, just add `--report-DNA-sequences`. 

### Exporting only specific gene clusters

#### Part 1: Choosing gene clusters by collection, bin, or name

You can export only the sequences for a specific <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> or <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> with the parameters `-C` or `-b` respectively. You also have the option to display the collections and bins available in your <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> with `--list-collections` or `--list-bins`

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;clusters &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                     &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
                                     &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/7.1/artifacts/genes&#45;fasta)</span> \
                                     &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> 
</div>

Alternatively, you can export the specific gene clusters by name, either by providing a single gene cluster ID or a file with one gene cluster ID per line. For example: 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;clusters &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                     &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
                                     &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/7.1/artifacts/genes&#45;fasta)</span> \
                                     &#45;&#45;gene&#45;cluster&#45;ids&#45;file gene_clusters.txt
</div>

where `gene_clusters.txt` contains the following:

    GC_00000618
    GC_00000643
    GC_00000729

#### Part 2: Choosing gene clusters by their attributes

These parameters are used to exclude gene clusters that don't reach certain thresholds and are applies on top of filters already applied (for example, you can use these to exclude clusters within a specific bin). 

Here is a list of the different filters that you can use to exclude some subsection of your gene clusters:

- min/max number of genomes that the gene cluster occurs in. 
- min/max number of genes from each genome. For example, you could exclude clusters that don't appear in every genome 3 times, or get single-copy genes by setting `max-num-genes-from-each-genome` to 1. 
- min/max [geometric homogenity index](http://merenlab.org/2016/11/08/pangenomics-v2/#geometric-homogeneity-index) 
- min/max [functional homogenity index](http://merenlab.org/2016/11/08/pangenomics-v2/#functional-homogeneity-index)
- min/max combined homogenity index 

For example, the following run on a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> that contains 50 genomes will report only the single-copy core genes with a functional homogenity index above 0.25:

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;clusters &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                     &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
                                     &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/7.1/artifacts/genes&#45;fasta)</span> \
                                     &#45;&#45;max&#45;num&#45;genes&#45;from&#45;each&#45;genome 1 \
                                     &#45;&#45;min&#45;num&#45;genomes&#45;gene&#45;cluster&#45;occurs 50 \
                                     &#45;&#45;min&#45;functional&#45;homogenity&#45;index 0.25 
</div>

You can also exclude genomes that are missing some number of the gene clusters that you're working with by using the paramter `--max-num-gene-clusters-missing-from-genome`. 

For each of these parameters, see the program's help menu for more information. 

### Fun with phylogenomics! 

To get a <span class="artifact-n">[concatenated-gene-alignment-fasta](/software/anvio/help/7.1/artifacts/concatenated-gene-alignment-fasta)</span> (which you can use to run <span class="artifact-n">[anvi-gen-phylogenomic-tree](/software/anvio/help/7.1/programs/anvi-gen-phylogenomic-tree)</span>), use the parameter `--concatenate-gene-clusters`

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;clusters &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                     &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
                                     &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/7.1/artifacts/genes&#45;fasta)</span> \
                                     &#45;&#45;concatenate&#45;gene&#45;clusters
</div>

Here, you also have the option to specify a specific aligner (or list the available aligners), as well as provide a NEXUS formatted partition file, if you so choose. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-sequences-for-gene-clusters.md) to update this information.


## Additional Resources


* [In action in the Anvi&#x27;o pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#scrutinizing-phylogenomics)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-sequences-for-gene-clusters) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
