---
layout: page
title: anvi-gen-gene-level-stats-databases [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-gen-gene-level-stats-databases
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program to compute genes databases for a ginen set of bins stored in an anvi&#x27;o collection. Genes databases store gene-level coverage and detection statistics, and they are usually computed and generated automatically when they are required (such as running anvi-interactive with `--gene-mode` flag). This program allows you to pre-compute them if you don&#x27;t want them to be done all at once.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[genes-db](../../artifacts/genes-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **generates a <span class="artifact-n">[genes-db](/software/anvio/help/main/artifacts/genes-db)</span>, which stores the coverage and detection values for all of the genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>.** 

This information is usually calculated when it's needed (for example when running <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> in genes mode), but this program lets you break this process into two steps. This way, you can easily change the parameters of <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> without having to recalculate the gene-level statistics. 

Given a <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> pair, as well as a <span class="artifact-n">[collection](/software/anvio/help/main/artifacts/collection)</span>, this program will calculate the stats for the genes in each of your <span class="artifact-n">[bin](/software/anvio/help/main/artifacts/bin)</span>s and give each bin its own <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> that includes this information. 

For example, if a <span class="artifact-n">[collection](/software/anvio/help/main/artifacts/collection)</span> called `GENE_COLLECTION` contained the bins `bin_0001`, `bin_0002`, and `bin_0003` and you ran:

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;gene&#45;level&#45;stats&#45;databases &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                                    &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                                    &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/main/artifacts/collection)</span> 
</div>

Then it will create a directory called `GENES` that contains three <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> called `GENE_COLLECTION-bin_0001.db`, `GENE_COLLECTION-bin_0002.db`, and `GENE_COLLECTION-bin_0003.db`. In terms of output, this program is similar to <span class="artifact-n">[anvi-split](/software/anvio/help/main/programs/anvi-split)</span>: each of these databases can now be treated as self-contained anvi'o projects but they also contain the gene-level information. Thus, you then could run <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> in genes mode on one of these profile databases. 

You also have the option to provide a list of <span class="artifact-n">[bin](/software/anvio/help/main/artifacts/bin)</span> (either as a file or as a string) to anlyze instead of a single <span class="artifact-n">[collection](/software/anvio/help/main/artifacts/collection)</span>. 

### Other Parameters

You can also change the definition of an outlier nucleotide position or switch calculations to use the [INSeq/Tn-Seq](https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/in-seq-tn-seq.html) statistical methods. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-gene-level-stats-databases.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-gene-level-stats-databases) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
