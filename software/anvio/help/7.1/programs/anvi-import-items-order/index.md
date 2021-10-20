---
layout: page
title: anvi-import-items-order [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-import-items-order
image:
  featurerelative: ../../../images/header.png
  display: true
---

Import a new items order into an anvi&#x27;o database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-items-order-txt](../../artifacts/misc-data-items-order-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[dendrogram](../../artifacts/dendrogram) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[phylogeny](../../artifacts/phylogeny) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[misc-data-items-order](../../artifacts/misc-data-items-order) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program, as one might think, allows you to import a <span class="artifact-n">[misc-data-items-order-txt](/software/anvio/help/7.1/artifacts/misc-data-items-order-txt)</span> to describe a specific order of items stored in a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>, <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>, or <span class="artifact-n">[genes-db](/software/anvio/help/7.1/artifacts/genes-db)</span>.

<div class="codeblock" markdown="1">
anvi&#45;import&#45;items&#45;order &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                        &#45;i <span class="artifact&#45;n">[misc&#45;data&#45;items&#45;order&#45;txt](/software/anvio/help/7.1/artifacts/misc&#45;data&#45;items&#45;order&#45;txt)</span>
</div>

It may also be nice to give it a good name, so that it's easy to find in the interface.

<div class="codeblock" markdown="1">
anvi&#45;import&#45;items&#45;order &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                        &#45;i <span class="artifact&#45;n">[misc&#45;data&#45;items&#45;order&#45;txt](/software/anvio/help/7.1/artifacts/misc&#45;data&#45;items&#45;order&#45;txt)</span> \
                        &#45;&#45;name ORDER_NAME
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-import-items-order.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-import-items-order) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
