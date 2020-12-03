---
layout: page
title: anvi-get-enriched-functions-per-pan-group [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program that takes a pangenome, and a categorical layers additional data item, and generates a table describing functions that are enriched in those groups. If requested, a functional occurrence table across genomes is also generated.

See **[program help menu](../../../vignette#anvi-get-enriched-functions-per-pan-group)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[functional-enrichment-txt](../../artifacts/functional-enrichment-txt)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[misc-data-layers](../../artifacts/misc-data-layers)</span> <span class="artifact-r">[pan-db](../../artifacts/pan-db)</span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db)</span></p>

## Usage


This program returns a **matrix of functions that are enriched within specific groups in your pangenome**. 

You provide a <span class="artifact-n">[pan-db](/software/anvio/help/main/artifacts/pan-db)</span> and <span class="artifact-n">[genomes-storage-db](/software/anvio/help/main/artifacts/genomes-storage-db)</span> pair, as well as a <span class="artifact-n">[misc-data-layers](/software/anvio/help/main/artifacts/misc-data-layers)</span> that stores categorical data, and the program will consider each of the categories their own 'pan-group'. It will then find functions that are enriched in that group (i.e., functions that are associated with gene clusters that are characteristic of the genomes in that group). It returns this output as a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/main/artifacts/functional-enrichment-txt)</span>

{:.notice}
Note that your <span class="artifact-n">[genomes-storage-db](/software/anvio/help/main/artifacts/genomes-storage-db)</span> must have at least one functional annotation source for this to work. 

This helps you highlight functions or pathways that separate a specific pan-group and determine the functional core of your pangenome. For example, in the *Prochlorococcus* pangenome (the one used in [the pangenomics tutorial, where you can find morei info about this program](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)), this program finds that `Exonuclease VII` is enriched in the low-light pan-group. The output file provides various statistics about how confident the program is in making this association. 

### How does it work? 

What this program does can be broken down into three steps: 

1. Determining the pan-groups. Firstly, the program uses a <span class="artifact-n">[misc-data-layers](/software/anvio/help/main/artifacts/misc-data-layers)</span> (containing categorical, not numerical, data) to split the pangenome into several groups. For example, in the pangenome tutorial, this was the low-light and high-light groups. 
2.  Determine the "functional associations" for each of your gene clusters. In short, this is collecting the functional annotations for all of the genes in each cluster and assigning the one that appears most frequently to represent the entire cluster. 
3. Looking at the functional associations and their relative levels of abundance across the pan-groups. Specifically, it looks at the level that a particular gene cluster's functional association is unique to a single pan-group and the percent of genomes it appears in in each pan-group. This is what is reported in the output matrix, a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/main/artifacts/functional-enrichment-txt)</span> 

If you're still curious, check out [Alon's behind the scenes post](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome), which goes into a lot more detail.

### Okay, cool. Let's talk parameters. 

Here is the simplest run of this program: 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;enriched&#45;functions&#45;per&#45;pan&#45;group &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/main/artifacts/pan&#45;db)</span>\
                                          &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/main/artifacts/genomes&#45;storage&#45;db)</span> \
                                          &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/main/artifacts/functional&#45;enrichment&#45;txt)</span> \ 
                                          &#45;&#45;category&#45;variable CATEGORY
</div>

The parameter `--category-variable` gives the name of the categorical <span class="artifact-n">[misc-data-layers](/software/anvio/help/main/artifacts/misc-data-layers)</span> that you want to use to define your pan-groups. Note that this will consider all items not in a category in their own 'ungrouped' pan-group; you can ignore those items with the flag `--exlcude-ungrouped` 

You can choose to not group together gene clusters with the same function by adding the parameter `--include-gc-identity-as-function` and setting the annotation source ot `IDENTITY`

#### Defining a Functional Annotation Source

You also have the option to only use functional annotations from a specific source, like so: 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;enriched&#45;functions&#45;per&#45;pan&#45;group &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/main/artifacts/pan&#45;db)</span>\
                                          &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/main/artifacts/genomes&#45;storage&#45;db)</span> \
                                          &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/main/artifacts/functional&#45;enrichment&#45;txt)</span> \ 
                                          &#45;&#45;category&#45;variable CATEGORY \ 
                                          &#45;&#45;annotation&#45;source COG_FUNCTION
</div>

Use the parameter `--list-annotation-sources` to list the available annotation sources in your <span class="artifact-n">[pan-db](/software/anvio/help/main/artifacts/pan-db)</span>. 

#### Additional Output 

You can also output a functional occurance table, which describes the number of times each of your functional associations occurs in each genome you're looking at. You can interact more with this data by using <span class="artifact-n">[anvi-matrix-to-newick](/software/anvio/help/main/programs/anvi-matrix-to-newick)</span>. 

You can find more information about this [here](http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-enriched-functions-per-pan-group.md) to update this information.


## Additional Resources


* [An example usage of the program in the context of the Prochlorococcus metapangenome from Delmont and Eren 2018 is included in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/)

* [An older version of this program was used in this publication by Benjamin Tully](https://www.nature.com/articles/s41467-018-07840-4)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-enriched-functions-per-pan-group) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
