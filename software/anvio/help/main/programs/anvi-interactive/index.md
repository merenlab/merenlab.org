---
layout: page
title: anvi-interactive [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Start an anvi&#39;o server for the interactive interface.

See **[program help menu](../../../vignette#anvi-interactive)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[collection](../../artifacts/collection)</span> <span class="artifact-p">[bin](../../artifacts/bin)</span> <span class="artifact-p">[interactive](../../artifacts/interactive)</span> <span class="artifact-p">[svg](../../artifacts/svg)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[single-profile-db](../../artifacts/single-profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[genes-db](../../artifacts/genes-db)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span> <span class="artifact-r">[view-data](../../artifacts/view-data)</span> <span class="artifact-r">[dendrogram](../../artifacts/dendrogram)</span> <span class="artifact-r">[phylogeny](../../artifacts/phylogeny)</span></p>

## Usage


Anvi-interactive opens the Anvi'o interactive interface, which is one of the most sophisticated parts of Anvi'o.

The most widely-known view gives you beautiful concentric circles of data, but the interface has many forms and vast functionality, from manual metagenomic binning (check out anvi-interactive in a metagenomic workflow [here](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-interactive)) to examining various stats about your data. However, in case you don't like circles, you can also display your data in a rectangle (as seen [here](http://merenlab.org/tutorials/interactive-interface/#lets-go-all-corners)).

In fact, the interface has many of its own blog posts, including a pretty comprehensive introductory tutorial [here](http://merenlab.org/tutorials/interactive-interface/) and a breakdown of its data types  [here](http://merenlab.org/2016/02/27/the-anvio-interactive-interface/).

Here, we'll go through *some* things that the Anvi'o interactive interface is capable of through this program. More information about most of this can be found by calling `anvi-interactive -h` or by checking out the additional resources at the bottom of this page.

## Running anvi-interactive on a profile database

One of the simplest ways to run the interactive interface (especially useful for manual binning) is just providing a profile database and contigs database:

<div class="codeblock" markdown="1">
anvi&#45;interactive &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span>
</div>

For the central tree to display correctly, you'll need to have run hierarchical clustering at some point while making your profile database (either during <span class="artifact-n">[anvi-merge](/software/anvio/help/main/programs/anvi-merge)</span>, or, if this is a <span class="artifact-n">[single-profile-db](/software/anvio/help/main/artifacts/single-profile-db)</span>, while running <span class="artifact-n">[anvi-profile](/software/anvio/help/main/programs/anvi-profile)</span>).

You'll get lovely concentric circles (or rectangles), each filled with data that was contained in your databases and that you are now free to interact with. See the page for the <span class="artifact-n">[interactive](/software/anvio/help/main/artifacts/interactive)</span> interface for more information.

### Running on a specific collection

You can also run anvi-interactive on a specific collection. When doing this, Anvi'o will calculate various information about each of your bins, and display the interface. Each item of your central plot will not represent a contig, but a bin within your collection.

<div class="codeblock" markdown="1">
anvi&#45;interactive &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                 &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/main/artifacts/collection)</span>
</div>

Since clustering is done here, you can also customize the linkage method and distance metric if desired.

See the note on this mode in [the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-interactive) for more information.

If instead you want to run the standard anvi'o interface, but only looking at contigs within a specific collection, use the tag `--collection-autoload`.

### Visualizing *gene coverages* instead of bins or contigs

You can also start the interface in "gene mode", in which each item of the central tree is a gene instead of a split or contig (or bin like in "collection mode").

A lot of the same functionality is available, including the inspect page and search options, however, you cannot store states or collections in this mode.

To initiate the visualization in gene mode you need the following:

<div class="codeblock" markdown="1">
anvi&#45;interactive &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                 &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/main/artifacts/collection)</span> \
                 &#45;b <span class="artifact&#45;n">[bin](/software/anvio/help/main/artifacts/bin)</span> \
                 &#45;&#45;gene&#45;mode
</div>

To see a practical example, please visit [this page](http://merenlab.org/tutorials/infant-gut/#the-gene-mode-studying-distribution-patterns-at-the-gene-level). For a real-world application of this feature, please take a look at this We previously have used this strategy to study gene coverages across metagenomes.

In this view you can order genes based on their distributions patterns across metagenomes (*without paying attention to their synteny*) or by ordering them based on their synteny in the genome (*without paying attention to their differential distribution*). [Figure 2 in this paper](https://peerj.com/articles/4320/) examples the latter, and [Figure 5 in this paper](https://stm.sciencemag.org/content/11/507/eaau9356) examples the former case, which is also shown below:

![](http://merenlab.org/images/gene-distribution-across-metagenomes.png)

## Running anvi-interactive in manual mode

You can initiate the anvi'o interactive interface in manual mode to run it on ad-hoc data (here is [a tutorial on this](http://merenlab.org/tutorials/interactive-interface/)).

Anvi'o interactive interface is initiated with the flag `--manual-mode` and then by providing *any* of the following types of files individually or together:

- a TAB-delimited tabular data,
- a NEWICK formatted tree,

When doing this kind of run, anvi'o does not expect you to have a profile database, but it still needs you to provide a name for it. Anvi'o will simply create an empty one for you so you can store your state or collections in it for reproducibility.

## Visualization Settings

In Anvi'o, the visualization settings at a given time are called a <span class="artifact-n">[state](/software/anvio/help/main/artifacts/state)</span>.

To open the interface in a specific state, you can use the `--state-autoload` flag or by importing a state using <span class="artifact-n">[anvi-import-state](/software/anvio/help/main/programs/anvi-import-state)</span>.

You can also customize various aspects of the interactive interface. For example, you can change the preselected view, title, and taxonomic level displayed (for example, showing the class name instead of the genus name). You can also hide outlier single nucleotide variations or open only a specific collection.

### Adding additional information

You can add any additional layers of your choice using the parameter `--additional-layers` and providing a file containing the information you want to appear as another layer. You could also choose to split non-single-copy gene HMM hits into their own layer with the `--split-hmm-layers` parameter.

If you want to add an entirely new view to the interface, you can do that too, as long as you provide a file containing all split names and their associated values. For more information, see the parameter `--additional-view`.

You can also provide the manual inputs even if you're using an Anvi'o database. For example, if you provide your own NEWICK formatted tree, you will have the option to display it instead of the one in your database.

## Other things

### Viewing your data

You can use this program to look at the available information in your databases, which is very convenient. For example, you can view all of the available

- views (using `--show-views`)
- states (using `--show-states`)
- collections (using `--list-collections`)

### Note for power users

You can also configure the server to your heart's content, set a password to make sure only those who has the password can work with the contents, skip initiating functions, and change any of the output paths, automatically save figures, etc.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-interactive.md) to update this information.


## Additional Resources


* [A beginners tutorial on anvi&#39;o interactive interface](http://merenlab.org/tutorials/interactive-interface/)

* [How to add more data to a display for layers and items](http://merenlab.org/2017/12/11/additional-data-tables/)

* [An overview of interactive data types](http://merenlab.org/2016/02/27/the-anvio-interactive-interface/)

* [Anvi&#39;o &#39;views&#39; demystified](http://merenlab.org/2017/05/08/anvio-views/)

* [Working with SVG files from the interactive interface](http://merenlab.org/2016/10/27/high-resolution-figures/)

* [Running remote anvi&#39;o interactive interfaces on your local computer](http://merenlab.org/2018/03/07/working-with-remote-interative/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-interactive) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
