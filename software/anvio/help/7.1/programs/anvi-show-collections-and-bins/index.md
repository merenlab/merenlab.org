---
layout: page
title: anvi-show-collections-and-bins [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-show-collections-and-bins
image:
  featurerelative: ../../../images/header.png
  display: true
---

A script to display collections stored in an anvi&#x27;o profile or pan database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


This program does not seem to provide any artifacts. Such programs usually print out some information for you to see or alter some anvi'o artifacts without producing any immediate outputs.


## Usage


This program tells you about the <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>s within a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> or <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>. 

Just run it like so 

<div class="codeblock" markdown="1">
anvi&#45;show&#45;collections&#45;and&#45;bins &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> 
</div>

and Anvi'o will output to your console the following information for each of the <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>s in the database: 

* The name and ID of the collection
* The number of <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>s within the collection, and each of their names
* The number of splits contained within those bins 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-show-collections-and-bins.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-show-collections-and-bins) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
