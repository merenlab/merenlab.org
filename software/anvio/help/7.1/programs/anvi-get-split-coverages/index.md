---
layout: page
title: anvi-get-split-coverages [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-get-split-coverages
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export splits and the coverage table from database.

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


<p style="text-align: left" markdown="1"><span class="artifact-p">[coverages-txt](../../artifacts/coverages-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program returns the nucleotide-level coverage data for a specific set of the splits or gene in your <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>. 

If you want to get the coverage data for all splits in your <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>, run <span class="artifact-n">[anvi-export-splits-and-coverages](/software/anvio/help/7.1/programs/anvi-export-splits-and-coverages)</span> with the flag `--splits-mode`. 

Simply provide a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> pair and specify which splits, or gene, you want to look at. You have three ways to do this: 

1.  Provide a single split name. (You can list all splits available with `--list-splits`)

<div class="codeblock" markdown="1">
anvi&#45;get&#45;split&#45;coverages &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                         &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                         &#45;o <span class="artifact&#45;n">[coverages&#45;txt](/software/anvio/help/7.1/artifacts/coverages&#45;txt)</span> \ 
                         &#45;&#45;split&#45;name Day17a_QCcontig9_split_00003
</div>


2. Provide both the name of a <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> and the <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> it is contained in. 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;split&#45;coverages &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                         &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                         &#45;o <span class="artifact&#45;n">[coverages&#45;txt](/software/anvio/help/7.1/artifacts/coverages&#45;txt)</span> \ 
                         &#45;b <span class="artifact&#45;n">[bin](/software/anvio/help/7.1/artifacts/bin)</span> \
                         &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>
</div>

You can list all collections available with `--list-collections` or all bins in a collection with `--list-bins`. Alternatively, you could run <span class="artifact-n">[anvi-show-collections-and-bins](/software/anvio/help/7.1/programs/anvi-show-collections-and-bins)</span> on your <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> to get a more comprehensive overview. 

3. Provide a gene caller id and a flanking size (bp).

<div class="codeblock" markdown="1">
anvi&#45;get&#45;split&#45;coverages &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                         &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                         &#45;o <span class="artifact&#45;n">[coverages&#45;txt](/software/anvio/help/7.1/artifacts/coverages&#45;txt)</span> \ 
                         &#45;&#45;gene&#45;caller&#45;id 25961 \
                         &#45;&#45;flank&#45;length 500
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-split-coverages.md) to update this information.


## Additional Resources


* [Using this program to generate split coverage visualizations across samples](http://merenlab.org/2019/11/25/visualizing-coverages/#visualize-only-the-coverage-of-a-split-across-samples)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-split-coverages) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
