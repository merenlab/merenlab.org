---
layout: page
title: anvi-compute-functional-enrichment [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This is a driver program for `anvi-script-enrichment-stats`, a script that computes enrichment scores and group associations for annotated entities (ie, functions, KEGG Modules) across groups of genomes or samples..

See **[program help menu](../../../../vignette#anvi-compute-functional-enrichment)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[functional-enrichment-txt](../../artifacts/functional-enrichment-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[kegg-metabolism](../../artifacts/kegg-metabolism) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[groups-txt](../../artifacts/groups-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-layers](../../artifacts/misc-data-layers) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program has multiple abilities. It can compute enriched functions across categories in a pangenome, enriched metabolic modules across groups of samples, or enriched functions across groups of genomes. To do this it relies on the script `anvi-script-enrichment-stats` by [Amy Willis](https://github.com/adw96).

Regardless of the situation, it returns a **matrix of things that are enriched within specific groups in your dataset**, as a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7/artifacts/functional-enrichment-txt)</span> file.

## Input option 1: Enriched functions in a pangenome

In this case, this program will return a matrix of functions that are enriched within specific groups in your pangenome.

You provide a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> and <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7/artifacts/genomes-storage-db)</span> pair, as well as a <span class="artifact-n">[misc-data-layers](/software/anvio/help/7/artifacts/misc-data-layers)</span> that stores categorical data, and the program will consider each of the categories their own 'pan-group'. It will then find functions that are enriched in that group (i.e., functions that are associated with gene clusters that are characteristic of the genomes in that group). It returns this output as a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7/artifacts/functional-enrichment-txt)</span>.

{:.notice}
Note that your <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7/artifacts/genomes-storage-db)</span> must have at least one functional annotation source for this to work.

This helps you highlight functions or pathways that distinguish a specific pan-group and determine the functional core of your pangenome. For example, in the *Prochlorococcus* pangenome (the one used in [the pangenomics tutorial, where you can find more info about this program](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)), this program finds that `Exonuclease VII` is enriched in the low-light pan-group. The output file provides various statistics about how confident the program is in making this association.

### How does it work?

What this program does can be broken down into three steps:

1. Determining the pan-groups. Firstly, the program uses a <span class="artifact-n">[misc-data-layers](/software/anvio/help/7/artifacts/misc-data-layers)</span> (containing categorical, not numerical, data) to split the pangenome into several groups. For example, in the pangenome tutorial, this was the low-light and high-light groups.
2.  Determine the "functional associations" for each of your gene clusters. In short, this is collecting the functional annotations for all of the genes in each cluster and assigning the one that appears most frequently to represent the entire cluster.
3. Looking at the functional associations and their relative levels of abundance across the pan-groups. Specifically, it looks at the level that a particular gene cluster's functional association is unique to a single pan-group and the percent of genomes it appears in in each pan-group. This is what is reported in the output matrix, a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7/artifacts/functional-enrichment-txt)</span>.

If you're still curious, check out [Alon's behind the scenes post](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome), which goes into a lot more detail.

### Basic usage

Here is the simplest run of this program:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span>\
                                   &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7/artifacts/genomes&#45;storage&#45;db)</span> \
                                   &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                                   &#45;&#45;category&#45;variable CATEGORY \
                                   &#45;&#45;annotation&#45;source FUNCTION_SOURCE
</div>

You must provide this program with a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> and its corresponding <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7/artifacts/genomes-storage-db)</span>. You must also provide an output file name.

The <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> must contain at least one categorical data layer in <span class="artifact-n">[misc-data-layers](/software/anvio/help/7/artifacts/misc-data-layers)</span>, and you must choose one of these categories to define your pan-groups with the `--category-variable` parameter. Note that by default any genomes not in a category will be ignored; you can instead include these in the analysis by using the flag `--include-ungrouped`.

The <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7/artifacts/genomes-storage-db)</span> must have at least one functional annotation source, and you must choose one of these sources with the `--annotation-source`. If you do not know which functional annotation sources are available in your <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7/artifacts/genomes-storage-db)</span>, you can use the `--list-annotation-sources` parameter to find out.

### Additional options

By default, gene clusters with the same functional annotation will be merged. But if you provide the `--include-gc-identity-as-function` parameter and set the annotation source to be 'IDENTITY', anvi'o will treat gene cluster names as functions and enable you to investigate enrichment of each gene cluster independently. This is how you do it:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span>\
                               &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7/artifacts/genomes&#45;storage&#45;db)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                               &#45;&#45;category&#45;variable CATEGORY \
                               &#45;&#45;annotation&#45;source IDENTITY \
                               &#45;&#45;include&#45;gc&#45;identity&#45;as&#45;function
</div>

To output a functional occurrence table, which describes the number of times each of your functional associations occurs in each genome you're looking at, use the `--functional-occurrence-table-output` parameter, like so:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span>\
                               &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7/artifacts/genomes&#45;storage&#45;db)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                               &#45;&#45;category&#45;variable CATEGORY \
                               &#45;&#45;annotation&#45;source FUNCTION_SOURCE \
                               &#45;&#45;functional&#45;occurrence&#45;table&#45;output FUNC_OCCURRENCE.TXT
</div>

You can interact more with this data file by using <span class="artifact-n">[anvi-matrix-to-newick](/software/anvio/help/7/programs/anvi-matrix-to-newick)</span>. Find more information about this output option [here](http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions).

## Input option 2: Enriched modules

This option computes enrichment scores for metabolic modules in groups of samples. In order to do this, you must already have estimated completeness of metabolic modules in your samples using <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> and obtained a "modules" mode output file (the default). You must provide that file to this program along with a <span class="artifact-n">[groups-txt](/software/anvio/help/7/artifacts/groups-txt)</span> file indicating which samples belong to which groups.

### How does it work?

1. Determining presence of modules. Each module in the "modules" mode output has a completeness score associated with it in each sample, and any module with a completeness score over a given threshold (set by `--module-completion-threshold`) will be considered to be present in that sample.
2. Examining the distribution of modules in each group of samples to compute an enrichment score for each module. This is done by fitting a generalized linear model (GLM) with a logit linkage function in `anvi-script-enrichment-stats`, and it produces a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7/artifacts/functional-enrichment-txt)</span> file.

### Basic usage

See <span class="artifact-n">[kegg-metabolism](/software/anvio/help/7/artifacts/kegg-metabolism)</span> for more information on the "modules" mode output format from <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span>, which you must provide with the `-M` flag. The sample names in this file must match those in the <span class="artifact-n">[groups-txt](/software/anvio/help/7/artifacts/groups-txt)</span> file, provided with `-G`. You must also provide the name of the output file.

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;M MODULES.TXT \
                               &#45;G <span class="artifact&#45;n">[groups&#45;txt](/software/anvio/help/7/artifacts/groups&#45;txt)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span>
</div>

### Additional parameters

The default completeness threshold for a module to be considered 'present' in a sample is 0.75 (75 percent). If you wish to change this, you can do so by providing a different threshold - as a number in the range (0, 1] - using the `--module-completion-threshold` parameter. For example:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;M MODULES.TXT \
                               &#45;G <span class="artifact&#45;n">[groups&#45;txt](/software/anvio/help/7/artifacts/groups&#45;txt)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                               &#45;&#45;module&#45;completion&#45;threshold 0.9
</div>

By default, the column containing sample names in your MODULES.TXT file will have the header `db_name`, but there are certain cases in which you might have them in a different column - for example, if you did not run <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> in multi-mode. In those cases, you can specify that a different column contains the sample names by providing its header with `--sample-header`. For example, if you sample names were in the `metagenome_name` column, you would do the following:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;M MODULES.TXT \
                               &#45;G <span class="artifact&#45;n">[groups&#45;txt](/software/anvio/help/7/artifacts/groups&#45;txt)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                               &#45;&#45;sample&#45;header metagenome_name
</div>

If you ran <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> on a bunch of extra samples but only want to include a subset of those samples in the <span class="artifact-n">[groups-txt](/software/anvio/help/7/artifacts/groups-txt)</span>, that is fine - by default any samples from the MODULES.TXT file that are missing from the <span class="artifact-n">[groups-txt](/software/anvio/help/7/artifacts/groups-txt)</span> will be ignored. However, there is also an option to include those missing samples in the analysis, as one big group called 'UNGROUPED'. To do this, you can use the --include-samples-missing-from-groups-txt parameter. Just be careful that if you are also using the --include-ungrouped flag (see below), any samples without a specified group in the <span class="artifact-n">[groups-txt](/software/anvio/help/7/artifacts/groups-txt)</span> will also be included in the 'UNGROUPED' group.

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;M MODULES.TXT \
                               &#45;G <span class="artifact&#45;n">[groups&#45;txt](/software/anvio/help/7/artifacts/groups&#45;txt)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                               &#45;&#45;include&#45;samples&#45;missing&#45;from&#45;groups&#45;txt
</div>


## Input option 3: Enriched functions in groups of genomes

You are not limited to computing functional enrichment in pangenomes, you can do it for regular genomes, too. This option takes either external or internal genomes (or both) which are organized into groups, and computes enrichment scores and associated groups for annotated functions in those genomes.

### How does it work?

This is similar to computing functional enrichment in pangenomes (as described above), but a bit simpler.

1. Counting functions. Gene calls in each genome are tallied according to their functional annotations from the given annotation source.
2. Looking at the functions and their relative levels of abundance across the groups of genomes. This again uses `anvi-script-enrichment-stats` to fit a GLM to determine A) the level that a particular functional annotation is unique to a single group and B) the percent of genomes it appears in in each group. This produces a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7/artifacts/functional-enrichment-txt)</span> file.

### Basic usage

You can provide either an <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span> file or an <span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span> file or both, but no matter what these files must contain a `group` column which indicates the group that each genome belongs to. Similar to option 1, you must also provide an annotation source from which to extract the functional annotations of interest. In the example below, we provide both types of input files.

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7/artifacts/internal&#45;genomes)</span>\
                               &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/7/artifacts/external&#45;genomes)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                               &#45;&#45;annotation&#45;source FUNCTION_SOURCE
</div>

### Additional Parameters

Also similar to option 1, you can get a tab-delimited matrix describing the occurrence (counts) of each function within each genome using the `--functional-occurrence-table-output` parameter:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7/artifacts/internal&#45;genomes)</span>\
                               &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/7/artifacts/external&#45;genomes)</span> \
                               &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7/artifacts/functional&#45;enrichment&#45;txt)</span> \
                               &#45;&#45;annotation&#45;source FUNCTION_SOURCE
                               &#45;&#45;functional&#45;occurrence&#45;table&#45;output FUNC_OCCURRENCE.TXT
</div>


## Parameters common to all options

If you provide the `--include-ungrouped` parameter, then genomes (or samples) without a group will be included from the analysis. (By default, these genomes/samples are ignored.) For the pangenome case, these genomes are those without a category in the provided `--category-variable`. For metabolic modules or the genomes in groups case, these samples/genomes are those with an empty value in the 'group' column (of either the <span class="artifact-n">[groups-txt](/software/anvio/help/7/artifacts/groups-txt)</span> or the <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span>/<span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span> files).


## More information on `anvi-script-enrichment-stats`

This program serves as the interface to `anvi-script-enrichment-stats`, an R script which performs an enrichment test on your input. You will find a brief description of how this script works in Alon's "Behind the Scenes" note in [the pangenomics tutorial](https://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome). Better yet, check out the methods section of Alon's paper, published in Genome Biology [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-compute-functional-enrichment.md) to update this information.


## Additional Resources


* [A description of the enrichment script run by this program can be found in Shaiber et al 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w)

* [An example of pangenome functional enrichment in the context of the Prochlorococcus metapangenome from Delmont and Eren 2018 is included in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-compute-functional-enrichment) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
