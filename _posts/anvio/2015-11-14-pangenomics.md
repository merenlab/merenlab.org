---
layout: post
title: "An anvi'o workflow for microbial pangenomics"
excerpt: "The user-friendly interface anvi'o provides to work with pangenomes."
modified: 2015-10-14
tags: []
categories: [anvio]
comments: true
authors: [tom, meren]
---

{% include _toc.html %}
Cultivation of closely related microorganisms, and the subsequent recovery of their genomes, revealed that not all isolates share the same functional traits for a given population. In other words, individual isolates do not fully echo the functional complexity of naturally occurring microbial populations, and reaching to a plateau of functional recovery is only possible through large genomic collections. Overlapping and differing functions among the genomes of closely related organisms led to the introduction of a new concept, "pangenomics" [(Tettelin et al., 2005)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1216834/), where genes across genomes are segregated into distinct groups: the core-genome, genes detected in multiple yet not all genomes, and finally isolate-specific genes. Pangenomic investigations are now widely used to dissect the functional traits of microorganisms, and to uncover environmentally- and clinically-important patterns of functioning.Although [multiple bioinformatics software](http://omictools.com/pangenomics-c1590-p1.html) are available to generate and/or visualize pangenomes, these solutions do not necessary offer flexible work environments, and hence limit the user's ability to interact with their data. 

That's why we are happy to demonstrate how [anvi'o](https://peerj.com/articles/1319/) can be used to process, visualize, and manipulate pangenomic data in its user-friendly environment. Although we are still [working](http://github.com/meren/anvio) on new modules for a fully automatized pangenomic workflow, the current anvi'o interface already provide its user with the opportunity to analyze pangenomes combined with a variety of contextual data, and to generate high-quality, publication-ready figures.

Our goal in this article is to describe original pangenomic analyses of publically available genomic data to demonstrate the anvi'o pangenomic workflow. Please don't hesitate to use the comments section for your suggestions, our [discussion group]({% post_url anvio/2015-10-14-anvio-discussion-group %}) for your technical questions. You can also keep an eye on our [public code repository](http://github.com/meren/anvio) for new releases, or to report issues.

Please note that you can find every file mentioned in this article, including run scripts to invoke anvi'o interactive interface here: [http://dx.doi.org/10.6084/m9.figshare.1603512](http://dx.doi.org/10.6084/m9.figshare.1603512).
## Step 1: Creating input files

You will find all files metnioned in this section in `01-ECOLI-PANGENOME/` directory of our figshare archive.

### Generating protein clusters for the data matrixFor our analysis of ten *Escherichia coli* (EC) genomes, we downloaded genbank files from [ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/](ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/).

{:.notice}
If you would like to repeat this analysis on your collection of genomes, you can use prodigal to identify open reading frames, and use the resulting genbank files instead.

We clustered amino acid sequences for 47,415 open reading frames in these genomes using [ITEP](https://price.systemsbiology.org/research/itep/) ("Integrated Toolkit for Exploration of Pan-genomes", a useful tool we were introduced by our colleague [Rika Anderson](https://twitter.com/RikaEAnderson)). ITEP identified a total of 7,795 protein clusters in our collection, where each cluster contained 1 to 102 proteins. 

The main purpose of this step is to create a basic matrix to connect genomes and protein clusters. Although we used ITEP, any other solution to cluster protein sequences could have been used to acquire this essential input file.

We stored this information as a matrix of protein clusters and their occurrences across 10 genomes (see `data.txt` in our figshare archive) which should look like the following table, where each cell shows the count of a given protein cluster (first column) in each genome (first row):

|contig|O157_H7-EDL933|CFT073|K-12_W3110|O157-H7-Sakai|SE11|SE15|BL21-DE3|K-12_MG1655|O103-H2-12009|O111-H-11128|
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1|18|8|0|18|0|0|0|0|23|35|
|2|19|9|0|18|0|0|0|0|23|33|
|3|2|1|6|1|1|4|56|6|1|1|
|4|10|2|2|11|2|1|3|1|10|12|
|5|9|3|2|11|2|1|1|2|10|10|
|6|0|0|7|1|1|3|29|7|1|2|
|7|5|6|4|4|7|5|3|4|6|5|
|8|14|6|0|13|0|5|0|0|6|3|
|9|9|3|0|10|4|1|1|0|9|9|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|7795|0|0|0|0|0|0|0|0|0|1|
### Creating a tree of protein clustersWe created a tree based on the distribution of protein clusters across genomes (`tree.newick`). This tree will be used as the centerpiece in the interactive interface you will see very soon:
{% highlight bash %}anvi-matrix-to-newick data.txt -o tree.newick{% endhighlight %}{:.notice}The creation of tree.newick will be quite challenging if your data.txt contains more than 30,000 protein clusters, in this case removing protein clusters that are detected only in one genome may help.### Creating a samples database
We created another tree from the same data:

{% highlight bash %}anvi-matrix-to-newick data.txt -o tree-for-samples-order.txt --transpose{% endhighlight %}

Note the `--transopse` flag. The resulting tree will describe the clustering of genomes with respect to protein clusters they harbor. We stored the resulting tree in `samples-order.txt` file, the format of which is explained [here](http://merenlab.org/2015/11/10/samples-db/#the-samples-order-file-format).

We also have another file that provides us with more information about our genomes. This file is called `samples-information.txt`, and it looks like this:

|samples|Nb_proteins|GC-content|
|:----|:----:|:----:|
|BL21-DE3|4192|50.8|
|CFT073|5364|50.5|
|O103-H2-12009|5049|50.7|
|O111-H-11128|4967|50.6|
|O157_H7-EDL933|5285|50.4|
|O157-H7-Sakai|5204|50.5|
|SE11|4673|50.8|
|SE15|4336|50.7|
|K-12_MG1655|4132|50.8|
|K-12_W3110|4213|50.8|

Using these two files, we then generate the samples database:

{% highlight bash %}anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o SAMPLES.db{% endhighlight %}

Please read [this blog post](http://merenlab.org/2015/11/10/samples-db/) to better make sense of anvi'o samples databases.
We are almost done!### Creating a mock FASTA fileIn the near future this step will not be necessary, but if you are reading this section, then you probably should do it.

We are going to create a mock FASTA file that contains entries for each protein to keep the interactive interface happy. For this, we simply take every protein cluster ID from our `data.txt` matrix, and add into a FASTA file with some mock sequences (clearly you can put in the actual sequences in this file, but the current version of anvi'o will not do much with them):

{% highlight bash %}
for i in `awk '{if(NR > 1) print $1}' data.txt`; do echo ">$i"; echo "ATCG"; done > fasta.fa{% endhighlight %}## Step 2: Visualizing the pangenome in the interactive interfaceAt this point we have our `data.txt` (the table for protein clusters), `tree.newick` (the tree for protein clusters), `SAMPLES.db` (anvi'o samples database to make sense of the contextual information about our genomes), and `fasta.fa` (our mock FASTA file that contains an entry for each protein cluster). It is time to visualize this pangenome, by invoking anvi'o interactive interface in manual mode:{% highlight bash %}anvi-interactive -f fasta.fa -p PROFILE.db -d data.txt -t tree.newick -s SAMPLES.db --state default_state --manual
{% endhighlight %}And voilà! Our simple pangenome of ten *E. coli* genomes are ready:<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/pan-genome-1.png"><img src="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/pan-genome-1.png" style="border: None;" /></a>
</div>

<div class="quotable">
Figure 1: Anvi'o pangenome display of ten *Escherichia coli* genomes with a total of 7,795 protein clusters. The inner tree (from `tree.newick`) organizes protein clusters based on their occurrence patterns across genomes. Layers around the tree display the presence/absence of protein clusters in a given genome (from `data.txt`). Note that the genome layers are organized based on the distribution of protein clusters (from `SAMPLES.db`). Finally, the number of proteins and GC-content (also from `SAMPLES.db`) is displayed to demonstrate the capabilities of anvi'o incorporating contextual data as a natural extention of the pangenome. </div>This figure demonstrates the pangenomic visualization capabilities of anvi'o in simplest terms. However, we didn't want to stop there, and also prepared a three-minute [video](https://youtu.be/RPjdl0bFCn0) to show how we generated this figure while demonstrating how the anvi'o interactive interface works (you should be able to read the text if you siwtch to full screen in HD):

<iframe width="560" height="315" src="https://www.youtube.com/embed/RPjdl0bFCn0" frameborder="0" allowfullscreen></iframe>
<div style="padding-bottom: 50px;"></div>
Thanks to the interactivity the interface provides, and the ability to extend the displayed information with additional simple matrices, we believe anvi'o fits quite nicely into the needs of analyses of pangenomes.## Combining multiple pangenomes for a multi-clade investigationHere we expand our analysis with more pangenomes, spanning beyond a single clade in the bacterial domain.In addition to the ten *Escherichia coli* genomes we analyzed in the previous section, we acquired genbank files for ten *Salmonella enterica*, and ten *Staphylococcus aureus* genomes. As you know, *E. coli* and *S. enterica* are in the same family, Enterobacteriaceae, and relatively distant from *S. aureus*.

Using ITEP, we organized the 119,511 genes identified in these genomes into 14,981 protein clusters, and generated all the relevant files as described in the first section. As an additional step, this time we also uploaded genomes into [RAST](http://rast.nmpdr.org/) to recover their metabolic profiles, and incorporated this contextual information into the metadata table before generating the `SAMPLES.db` file. All relevant files for this analysis are available in the directory `02-ECOLI-SAUREUS-SENTERICA-PANGENOME` in our figshare archive.In the interactive interface we colored sample layers based on their taxonomic affiliations, and added some extra space between groups to increase the visual clarity. We also selected five groups of protein clusters that represent important pangenomic features, and annotated these selections with numbers in the figure:<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/pan-genome-2.png"><img src="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/pan-genome-2.png" style="border: None;" /></a>
</div><div class="quotable">Figure 2: Pangenomic visualization of 14,981 protein clusters identified in ten <em>Escherichia coli</em>, ten <em>Salmonella enterica</em>, and ten <em>Staphylococcus aureus</em> genomes with anvi'o. The taxonomical affiliation, GC-content, number of proteins, and metabolic profiles of each genome is shown on the side. Finally, the tree contains selections of protein cluster groups displaying pangenomically relevant patterns: The core-pangenome, proteins characteristic to <em>E. Coli</em>, <em>S. enterica</em>, family Enterobacteriaceae, and <em>S. aureus</em>.
</div>As expected, the genomic clustering is consistent with their taxonomical affiliations. *E. coli* and *S. enterica* genomes share a large portion of protein clusters. On the other hand, only a small number of protein clusters are common to all 30 genomes. While some interesting traits (e.g., occurrence of phages) could be further explored, this is beyond the scope of our little blog post.## Final remarksApart from the protein clustering step, we completed these two pangenomic analyses within hours, including selecting core and peripheral genes, and creating figures in this post. In fact, writing this blog post took more time than both steps combined. Our analyses can be fully reproduced by downloading related anvi'o files from [here](http://dx.doi.org/10.6084/m9.figshare.1603512). We did this analyses using anvi'o version 1.2.1. You can have interactive access to these figures in a matter of minutes by running anvi'o through our [docker image]({% post_url anvio/2015-08-22-docker-image-for-anvio %}).

## Quick workflow for experimented users

The following script describes an overview of steps detailed in this post:

<script src="https://gist.github.com/meren/e99c46484ca5a479cf0c.js"></script>

Please don't hesitate to use the comments section for your suggestions, our [discussion group]({% post_url anvio/2015-10-14-anvio-discussion-group %}) for your technical questions. You can keep an eye on our [public code repository](http://github.com/meren/anvio) for new releases.