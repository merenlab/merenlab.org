---
layout: page
title: anvi-export-items-order [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-export-items-order
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export an item order from an anvi&#x27;o database.

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


<p style="text-align: left" markdown="1"><span class="artifact-p">[misc-data-items-order-txt](../../artifacts/misc-data-items-order-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[dendrogram](../../artifacts/dendrogram) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[phylogeny](../../artifacts/phylogeny) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span></p>


## Usage


This program, as one might think, allows you to export a <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7.1/artifacts/misc-data-items-order)</span>, outputing a <span class="artifact-n">[misc-data-items-order-txt](/software/anvio/help/7.1/artifacts/misc-data-items-order-txt)</span>. 

You can export one of the item orders in a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> or <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> as follows: 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;items&#45;order &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                        &#45;&#45;name cov
</div>

The `cov` here refers to the tree that is generated using only differential coverage. Almost all anvi'o profile databases will also have available an items-order based on the tetranucleotide frequency called `tnf`, and one based on both called `tnf-cov`. 

However, to list the item orders available in this database, just don't include the name flag.  

<div class="codeblock" markdown="1">
anvi&#45;export&#45;items&#45;order &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> 
</div>

You'll get a `Config Error` that will tell you what item orders are available. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-items-order.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-items-order) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
