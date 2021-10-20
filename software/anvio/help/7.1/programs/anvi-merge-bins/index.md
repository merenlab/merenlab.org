---
layout: page
title: anvi-merge-bins [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-merge-bins
image:
  featurerelative: ../../../images/header.png
  display: true
---

Merge a given set of bins in an anvi&#x27;o collection.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span></p>


## Can provide


This program does not seem to provide any artifacts. Such programs usually print out some information for you to see or alter some anvi'o artifacts without producing any immediate outputs.


## Usage


This program **merges two or more <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>s together** into a single <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>.

To run this program, the bins that you want to merge must be contained within a single <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>. Just provide the collection name, the <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> or <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> you're working with, the bins that you want to merge, and the name of the output bin. 

To check what collections and bins are contained in a database, you can either run this program with the flag `--list-collections` or `--list-bins`, or you can run <span class="artifact-n">[anvi-show-collections-and-bins](/software/anvio/help/7.1/programs/anvi-show-collections-and-bins)</span>.

For example, if you wanted to merge the bins `first_third`, `middle_third`, and `last_third` in a pan-db into a single bin called `complete_bin`, just run 

<div class="codeblock" markdown="1">
anvi&#45;merge&#45;bins &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
                &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> \
                &#45;b first_third,middle_third,last_third \
                &#45;B complete_bin
</div>

Now your collection will contain the bin `complete_bin` and the original bins will be gone forever (unless you had run<span class="artifact-n">[anvi-summarize](/software/anvio/help/7.1/programs/anvi-summarize)</span>, <span class="artifact-n">[anvi-export-collection](/software/anvio/help/7.1/programs/anvi-export-collection)</span>, or a similar program beforehand)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-merge-bins.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-merge-bins) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
