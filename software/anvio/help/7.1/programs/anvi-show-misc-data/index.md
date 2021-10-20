---
layout: page
title: anvi-show-misc-data [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-show-misc-data
image:
  featurerelative: ../../../images/header.png
  display: true
---

Show all misc data keys in all misc data tables.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


This program does not seem to provide any artifacts. Such programs usually print out some information for you to see or alter some anvi'o artifacts without producing any immediate outputs.


## Usage


This program **lists the additional data** that is stored within a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>, <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> or <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. This is data that can be imported with <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span> and is displayed in the interactive interface. 

When run, this program will output to the terminal a list of all additional data tables that are stored within the database. If you want to export a specific element of these as a text file, see <span class="artifact-n">[anvi-export-misc-data](/software/anvio/help/7.1/programs/anvi-export-misc-data)</span>. 

### What is displayed? 

When running on a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> or <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>, the output will display the following types of data:

- <span class="artifact-n">[misc-data-items](/software/anvio/help/7.1/artifacts/misc-data-items)</span> 
- <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span>
- <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7.1/artifacts/misc-data-layer-orders)</span> (by default, this will include orders like `abundance` and `mean_coverage (newick)`)

When running on a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, the output will display the following types of data:

- <span class="artifact-n">[misc-data-nucleotides](/software/anvio/help/7.1/artifacts/misc-data-nucleotides)</span> 
- <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/7.1/artifacts/misc-data-amino-acids)</span> 

These have no default values and will only contain data that has been imported with <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>. 

You also have the option to specify a specific kind of additional data table with `-t`. For example, to view only <span class="artifact-n">[misc-data-items](/software/anvio/help/7.1/artifacts/misc-data-items)</span> in a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>, just call

<div class="codeblock" markdown="1">
anvi&#45;show&#45;misc&#45;data &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                    &#45;t items 
</div>

Similarly to importing and exporting additional data tables, you can also focus on a specific data group with the parameter `-D`.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-show-misc-data.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-show-misc-data) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
