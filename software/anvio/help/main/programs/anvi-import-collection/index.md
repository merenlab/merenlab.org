---
layout: page
title: anvi-import-collection [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-import-collection
image:
  featurerelative: ../../../images/header.png
  display: true
---

Import an external binning result into anvi&#x27;o.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection-txt](../../artifacts/collection-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>


## Usage


This program, as one might think, allows you to import a <span class="artifact-n">[collection](/software/anvio/help/main/artifacts/collection)</span>. This allows you to easily import any binning that you've already done into a <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span>, since the <span class="artifact-n">[bin](/software/anvio/help/main/artifacts/bin)</span>s within that collection will be carried over. 

This information (in the form of a <span class="artifact-n">[collection-txt](/software/anvio/help/main/artifacts/collection-txt)</span>) can either come from another Anvi'o project (using <span class="artifact-n">[anvi-export-collection](/software/anvio/help/main/programs/anvi-export-collection)</span>) or you can get the coverage and sequence composion of your data using <span class="artifact-n">[anvi-export-splits-and-coverages](/software/anvio/help/main/programs/anvi-export-splits-and-coverages)</span> to bin your contigs with software outside of Anvi'o, then import that data into your database with this program. 

You can run this program like so: 

<div class="codeblock" markdown="1">
anvi&#45;import&#45;collection &#45;C my_bins.txt \
                        &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                        &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> 
</div>

This will import the collection indicated in `my_bins.txt` into your <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span>. 

`my_bins.txt` should be a tab-delimited file where the first column lists a split name and the second lists the bin that it is placed in. You can see an example of this [here](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_splits.txt). 

You can also provide this information by listing your contigs instead of your splits (like [this](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_contigs.txt)). Just add the `--contigs-mode` tag. 

You can also provide an information file to describe the source and/or colors of your bins. [This file](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/example_bins_info_file.txt) is an example of such an information file. 




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-import-collection.md) to update this information.


## Additional Resources


* [Another description as part of the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-collection)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-import-collection) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
