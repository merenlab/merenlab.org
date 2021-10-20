---
layout: page
title: anvi-pan-genome [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-pan-genome
image:
  featurerelative: ../../../images/header.png
  display: true
---

An anvi&#x27;o program to compute a pangenome from an anvi&#x27;o genome storage.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-items-order](../../artifacts/misc-data-items-order) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **finds, clusters, and organizes the genes** within a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> to create a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>. 

This is the program that does the brunt of the work when running a pangenomic workflow. Check out [the pangenomic tutorial](http://merenlab.org/2016/11/08/pangenomics-v2) for a more in-depth overview of the contents of this page and the capabilities of a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>. 

### Before running this program

Before running this program, you'll want to make sure your dependencies are all set, since this program requires some aditional dependencies than base anvi'o. If the following command runs without errors, then you're all good. 

<div class="codeblock" markdown="1">
anvi&#45;self&#45;test &#45;&#45;suite pangenomics
</div>

If that command doesn't run smoothly, check out [this page](http://merenlab.org/2016/11/08/pangenomics-v2/#dependencies).

### What this program does

This program finds and organizes your gene clusters to give you all of the data that is displayed when you run <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/7.1/programs/anvi-pan-genome)</span>. Almost all of the work described in [this gif that explains the common steps involved in pangenomics](http://merenlab.org/momics/#pangenomics) is done by this program. 

In a little more detail, this program will do three major things for you:

* Calculate the similarity between the all of the gene calls in all of the genomes in your <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span>. By default this uses [DIAMOND](https://www.wsi.uni-tuebingen.de/lehrstuehle/algorithms-in-bioinformatics/software/diamond/) to do this, but Meren strongly recommends that you use the `--use-ncbi-blast` flag to use [`blastp`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) instead.  

    *   When doing this, this will look at every genome in your <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> (unless you use `--genome-names`) and will use every gene call, whether or not they are complete (unless you used `--exclude-partial-gene-calls`).   
    
    *   After doing this, it will use the minbit heuristic (originally from ITEP ([Benedict et al., 2014](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-8)) to throw out weak matches. This removes a lot of noise before clustering. 
    
* Use the [MCL](http://micans.org/mcl/) algorithm to identify clusters in your search results.  

* Organize your gene clusters and genomes using their `euclidean` distance and `ward` linkage. 

This program is very smart, and if you're already run it, it will try to use the data that it's already calculated. This way you can change smaller parameters without all of the run time. However, this also means you need to tell it to rerun the process (if that's what you want) with the flag `--overwrite-output-destinations`. 

### Cool. How about some examples and specific parameters?

Who doesn't love a good example? The simplest way to run this is as follows:

<div class="codeblock" markdown="1">
anvi&#45;pan&#45;genome &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span>
</div>

But there are many parameters you can alter to your liking. For example, here's a run that specifies that it wants to use NCBI's blastp to find sequence similarities and muscle to align genes and defines its output 

<div class="codeblock" markdown="1">
anvi&#45;pan&#45;genomes &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                 &#45;&#45;align&#45;with muscle \
                 &#45;&#45;use&#45;ncbi&#45;blast \ 
                 &#45;n MY_PROJECT_NAME \
                 &#45;&#45;description description.txt \
                 &#45;o PATH/TO/<span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> 
</div>

Here's another example that only looks at the complete gene calls within a subset of the genomes, eliminates gene clusters that only have hits in a single genome, and uses DIAMOND but with the sensitive setting enabled:

<div class="codeblock" markdown="1">
anvi&#45;pan&#45;genomes &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                 &#45;n MY_PROJECT_NAME \
                 &#45;&#45;genome&#45;names GENOME_1,GENOME_2,GENOME_3 \
                 &#45;&#45;exclude&#45;partial&#45;gene&#45;calls \ 
                 &#45;&#45;min&#45;occurance 2 \
                 &#45;&#45;sensitive \
                 &#45;o PATH/TO/<span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> 
</div>

Some other parameters available to you allow you to  

- Change the minimum minbit value to elimate weak matches. The default value is 0.5.
- Change the MCL inflation parameter. The default value is 2. 
- Specify a minimum percent identity between two sequences to give that link an edge in MCL clustering. 
- Skip or speed up the calculation of homogeneity values for your clusters
- Enable multithreading with `-T`


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-pan-genome.md) to update this information.


## Additional Resources


* [A tutorial on pangenomics](http://merenlab.org/2016/11/08/pangenomics-v2/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-pan-genome) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
