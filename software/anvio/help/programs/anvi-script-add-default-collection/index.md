---
layout: page
title: anvi-script-add-default-collection [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

A script to add a &#39;DEFAULT&#39; collection in an anvi&#39;o pan or profile database with a bin named &#39;EVERYTHING&#39; that describes all items available in the profile database.

See **[program help menu](../../../vignette#anvi-script-add-default-collection)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[collection](../../artifacts/collection)</span> <span class="artifact-p">[bin](../../artifacts/bin)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


This program adds a new <span class="artifact-n">[collection](/software/anvio/help/artifacts/collection)</span> and <span class="artifact-n">[bin](/software/anvio/help/artifacts/bin)</span> to your <span class="artifact-n">[pan-db](/software/anvio/help/artifacts/pan-db)</span> or <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> pair. This collection and bin will both contain all of your contigs. 

This way, you can perform collection and bin specfic operations without having to bin anything yourself. For example, running <span class="artifact-n">[anvi-interactive](/software/anvio/help/programs/anvi-interactive)</span> in gene-mode requires you to specify a collection and bin (as is done [in the Infant Gut Tutorial](http://merenlab.org/tutorials/infant-gut/#the-gene-mode-studying-distribution-patterns-at-the-gene-level)). 

By default, the collection is named `DEFAULT` and the bin is named `EWVERYTHING`, but you can change these names with the `-C` and `-b` parameters respectively. 

Here is an example run on a <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span>: 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;add&#45;default&#45;collection &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \ 
                                   &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/artifacts/profile&#45;db)</span> \ 
                                   &#45;C MY_COLLECTION \
                                   &#45;b MY_BIN 
</div>

Once this is run, your profile database will contain a collection called `MY_COLLECTION` with a single bin (called `MY_BIN`) which contains all of your contigs. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-add-default-collection.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-add-default-collection) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
