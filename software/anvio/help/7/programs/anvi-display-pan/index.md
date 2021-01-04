---
layout: page
title: anvi-display-pan [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Start an anvi&#x27;o server to display a pan-genome.

See **[program help menu](../../../../vignette#anvi-display-pan)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[svg](../../artifacts/svg) <img src="../../images/icons/SVG.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program **displays the contents of a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> in the [anvi'o interactive interface](http://merenlab.org/2016/02/27/the-anvio-interactive-interface//#using-the-anvio-interactive-interface), much like <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span>.**

Like you can see in the [pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#displaying-the-pan-genome), this opens a window of the interactive interface where each item is a gene cluster and each layer represents one of your genomes. 

### A general run 

You can run it with only two parameters: 

<div class="codeblock" markdown="1">
anvi&#45;display&#45;pan &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span> \
                 &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7/artifacts/genomes&#45;storage&#45;db)</span> 
</div>

There are several default layer orders to choose from, including organizing based on gene cluster presence/absense or gene cluster frequency. These will both group your core gene clusters and singletons separately. 

Beyond that, there are many different settings you can change in the side panel of the interface and you can import various additional data (primarily with the program <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7/programs/anvi-import-misc-data)</span>). Once you're happy with the data displayed in the interface (and the prettiness of that data), you can  save those preferences in a <span class="artifact-n">[state](/software/anvio/help/7/artifacts/state)</span>. 

### I want MORE data displayed 

There are several other data types you can additionally choose to display in this program. Namely, you can add:

- a title (very fancy I know) using `--title` 
- a NEWICK formatted tree (or import it as a <span class="artifact-n">[misc-data-items-order-txt](/software/anvio/help/7/artifacts/misc-data-items-order-txt)</span> with <span class="artifact-n">[anvi-import-items-order](/software/anvio/help/7/programs/anvi-import-items-order)</span> or as a <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7/artifacts/misc-data-layer-orders)</span> with <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7/programs/anvi-import-misc-data)</span>). 
- view data in a tab-delimited file
- an additional view (provide this in a tab-delimited matrix where each column corresponds to a sample and each row corresponds to a gene cluster)
- an additional layer in the form of a <span class="artifact-n">[misc-data-layers-txt](/software/anvio/help/7/artifacts/misc-data-layers-txt)</span> (or import it into your <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> with <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7/programs/anvi-import-misc-data)</span>

### How to minimize mouse clicks 

Wondering how to autoload specific aspects of the interface? You're in the right place. 

You have the option to specify quite a few aspects of the interface through the parameters to save you those sweet mouse clicks. 

- You can specify which view to start the interface with. Check which views are available with `--list-views`. 
- You can load a specific <span class="artifact-n">[state](/software/anvio/help/7/artifacts/state)</span> (either a previous state or a state that you've imported with <span class="artifact-n">[anvi-import-state](/software/anvio/help/7/programs/anvi-import-state)</span>). Check which states are available with the flag `--list-states`. 
- You can also load a specific <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> with `--collection-autoload`. To check which collections are availible, use `--list-collections`. 

### Other parameters 

You can also skip processes like intializing functions or automatically ordering your items to save time, as well as configure the server to your heart's content. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-display-pan.md) to update this information.


## Additional Resources


* [See this program in action on the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#displaying-the-pan-genome)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-display-pan) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
