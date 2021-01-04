---
layout: page
title: anvi-import-collection [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Import an external binning result into anvi&#x27;o.

See **[program help menu](../../../../vignette#anvi-import-collection)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection-txt](../../artifacts/collection-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program, as one might think, allows you to import a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>. This allows you to easily import any binning that you've already done into a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>, since the <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span>s within that collection will be carried over. 

This information (in the form of a <span class="artifact-n">[collection-txt](/software/anvio/help/7/artifacts/collection-txt)</span>) can either come from another Anvi'o project (using <span class="artifact-n">[anvi-export-collection](/software/anvio/help/7/programs/anvi-export-collection)</span>) or you can get the coverage and sequence composion of your data using <span class="artifact-n">[anvi-export-splits-and-coverages](/software/anvio/help/7/programs/anvi-export-splits-and-coverages)</span> to bin your contigs with software outside of Anvi'o, then import that data into your database with this program. 

You can run this program like so: 

<div class="codeblock" markdown="1">
anvi&#45;import&#45;collection &#45;C my_bins.txt \
                        &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                        &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> 
</div>

This will import the collection indicated in `my_bins.txt` into your <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. 

`my_bins.txt` should be a tab-delimited file where the first column lists a split name and the second lists the bin that it is placed in. You can see an example of this [here](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_splits.txt). 

You can also provide this information by listing your contigs instead of your splits (like [this](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_contigs.txt)). Just add the `--contigs-mode` tag. 

You can also provide an information file to describe the source and/or colors of your bins. [This file](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/example_bins_info_file.txt) is an example of such an information file. 




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-import-collection.md) to update this information.


## Additional Resources


* [Another description as part of the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-collection)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-import-collection) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
