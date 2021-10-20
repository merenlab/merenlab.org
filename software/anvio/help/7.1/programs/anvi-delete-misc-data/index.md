---
layout: page
title: anvi-delete-misc-data [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-delete-misc-data
image:
  featurerelative: ../../../images/header.png
  display: true
---

Remove stuff from &#x27;additional data&#x27; or &#x27;order&#x27; tables for either items or layers in either pan or profile databases. OR, remove stuff from the &#x27;additional data&#x27; tables for nucleotides or amino acids in contigs databases.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-items](../../artifacts/misc-data-items) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-layers](../../artifacts/misc-data-layers) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-layer-orders](../../artifacts/misc-data-layer-orders) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-nucleotides](../../artifacts/misc-data-nucleotides) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-amino-acids](../../artifacts/misc-data-amino-acids) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Can provide


This program does not seem to provide any artifacts. Such programs usually print out some information for you to see or alter some anvi'o artifacts without producing any immediate outputs.


## Usage


After you've added misc-data of some kind (<span class="artifact-n">[misc-data-items](/software/anvio/help/7.1/artifacts/misc-data-items)</span>, <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span>, <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7.1/artifacts/misc-data-layer-orders)</span>, <span class="artifact-n">[misc-data-nucleotides](/software/anvio/help/7.1/artifacts/misc-data-nucleotides)</span>, or <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/7.1/artifacts/misc-data-amino-acids)</span>) using <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>, you can **delete that data and remove it from the interactive interface** using this program. 

This program will release your data into the ether, never to be seen again. If you would like to first export it into a text file (so that it can be seen again), you can do so with <span class="artifact-n">[anvi-export-misc-data](/software/anvio/help/7.1/programs/anvi-export-misc-data)</span>. 

This program only works on data that is listed as an available key (most often because it was previously imported by the user). To view available keys, call either

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;misc&#45;data &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                      &#45;&#45;target&#45;data&#45;table items|layers|layer_orders \
                      &#45;&#45;list&#45;available&#45;keys
</div>

or 

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;misc&#45;data &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;target&#45;data&#45;table nucleotides|amino_acids \
                      &#45;&#45;list&#45;available&#45;keys
</div>

where you choose the appropriate option for the `taget-data-table`. 

If your misc-data is associated with a specific data group, you can provide that data group to this program with the `-D` flag. 

## Data types you can delete 

### From a pan-db or profile-db: items, layers, layer orders

**From a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> or <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>, you can delete**

- items data (<span class="artifact-n">[misc-data-items](/software/anvio/help/7.1/artifacts/misc-data-items)</span>) 

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;misc&#45;data &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                      &#45;&#45;target&#45;data&#45;table items \
                      &#45;&#45;keys&#45;to&#45;remove key_1
</div>

- layers data (<span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span>)

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;misc&#45;data &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
                      &#45;&#45;target&#45;data&#45;table layers \
                      &#45;&#45;keys&#45;to&#45;remove key_2,key_3
</div>

- layer orders data (<span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7.1/artifacts/misc-data-layer-orders)</span>)

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;misc&#45;data &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                      &#45;&#45;target&#45;data&#45;table layer_orders \
                      &#45;&#45;keys&#45;to&#45;remove key_4
</div>

### From a contigs-db: nucleotide and amino acid data

**From a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, you can delete**

- nucleotide data (<span class="artifact-n">[misc-data-nucleotides](/software/anvio/help/7.1/artifacts/misc-data-nucleotides)</span>)

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;misc&#45;data &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;target&#45;data&#45;table nucleotides \
                      &#45;&#45;keys&#45;to&#45;remove key_1
</div>

- amino acid data (<span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/7.1/artifacts/misc-data-amino-acids)</span>)

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;misc&#45;data &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;target&#45;data&#45;table amino_acids \
                      &#45;&#45;keys&#45;to&#45;remove key_2
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-delete-misc-data.md) to update this information.


## Additional Resources


* [Working with anvi&#x27;o additional data tables](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-delete-misc-data) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
