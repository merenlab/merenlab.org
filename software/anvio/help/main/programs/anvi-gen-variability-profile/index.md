---
layout: page
title: anvi-gen-variability-profile [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate a table that comprehensively summarizes the variability of nucleotide, codon, or amino acid positions. We call these single nucleotide variants (SNVs), single codon variants (SCVs), and single amino acid variants (SAAVs), respectively.

See **[program help menu](../../../vignette#anvi-gen-variability-profile)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[variability-profile-txt](../../artifacts/variability-profile-txt)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[structure-db](../../artifacts/structure-db)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span> <span class="artifact-r">[variability-profile](../../artifacts/variability-profile)</span> <span class="artifact-r">[splits-txt](../../artifacts/splits-txt)</span></p>

## Usage


This program takes the variability data stored within a <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> and compiles it from across samples into a single matrix that comprehensively describes your SNVs, SCVs or SAAVs (a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/main/artifacts/variability-profile-txt)</span>).

By default, this program will just put the variability data into your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>. If you would like a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/main/artifacts/variability-profile-txt)</span> file, just provide the `-o` parameter with a filename. 

This program is described on [this blog post](http://merenlab.org/2015/07/20/analyzing-variability/#the-anvio-way), so take a look at that for more details. 

## Let's talk parameters 

Here is a basic run with no bells or whisles: 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                             &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \ 
                             &#45;C DEFAULT \
                             &#45;b EVERYTHING
</div>

Note that this program requires you to specify a subset of the databases that you want to focus on, so to focus on everything in the databases, run <span class="artifact-n">[anvi-script-add-default-collection](/software/anvio/help/main/programs/anvi-script-add-default-collection)</span> and use the resulting <span class="artifact-n">[collection](/software/anvio/help/main/artifacts/collection)</span> and <span class="artifact-n">[bin](/software/anvio/help/main/artifacts/bin)</span>, as shown above. 

You can add structural annotations by providing a <span class="artifact-n">[structure-db](/software/anvio/help/main/artifacts/structure-db)</span>. 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                             &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                             &#45;C DEFAULT \
                             &#45;b EVERYTHING \
                             &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/main/artifacts/structure&#45;db)</span> 
</div>

### Focusing on a subset of the input 

Instead of focusing on everything (providing the collection `DEFAULT` and the bin `EVERYTHING`), there are three ways to focus on a subset of the input: 

1. Provide a list of gene caller IDs (as a parameter with the flag `--gene-caller-ids` as shown below, or as a file with the flag `--genes-of-interest`)

    <div class="codeblock" markdown="1">
    anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                                 &#45;&#45;gene&#45;caller&#45;ids GENE_1,GENE_2,GENE_3
    </div>

2. Provide a <span class="artifact-n">[splits-txt](/software/anvio/help/main/artifacts/splits-txt)</span> to focus only on a specific set of splits. 

    <div class="codeblock" markdown="1">
    anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                                 &#45;&#45;splits&#45;of&#45;intest <span class="artifact&#45;n">[splits&#45;txt](/software/anvio/help/main/artifacts/splits&#45;txt)</span>
    </div>
    
3. Provide some other <span class="artifact-n">[collection](/software/anvio/help/main/artifacts/collection)</span> and <span class="artifact-n">[bin](/software/anvio/help/main/artifacts/bin)</span>. 

    <div class="codeblock" markdown="1">
    anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \ 
                                 &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/main/artifacts/collection)</span> \
                                 &#45;b <span class="artifact&#45;n">[bin](/software/anvio/help/main/artifacts/bin)</span>
    </div>

### Additional Ways to Focus the Input 

When providing a <span class="artifact-n">[structure-db](/software/anvio/help/main/artifacts/structure-db)</span>, you can also limit your analysis to only genes that have structures in your database. 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                             &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                             &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/main/artifacts/structure&#45;db)</span> \
                             &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/main/artifacts/collection)</span> \
                             &#45;&#45;only&#45;if&#45;structure
</div>

You can also choose to look at only data from specific samples by providing a file with one sample name per line. For example

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                             &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                             &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/main/artifacts/collection)</span> \
                             &#45;&#45;samples&#45;of&#45;interest my_samples.txt
</div>

where `my_samples.txt` looks like this:

    DAY_17A
    DAY_18A
    DAY_22A
    
### SNVs vs. SCVs vs. SAAVs 

Which one you're analyzing depends entirely on the `engine` parameter, which you can set to `NT` (nucleotides), `CDN` (codons), or `AA` (amino acids). The default value is nucleotides. Note that to analyze SCVs or SAAVs, you'll have needed to use the flag `--profile-SCVs` when you ran <span class="artifact-n">[anvi-profile](/software/anvio/help/main/programs/anvi-profile)</span> or <span class="artifact-n">[anvi-merge](/software/anvio/help/main/programs/anvi-merge)</span>. 

For example, to analyze SAAVs, run 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                             &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                             &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/main/artifacts/structure&#45;db)</span> \
                             &#45;&#45;engine AA
</div>

When analyzing single codon variants, you can choose to skip computing synonymity to save on run time, as so: 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;variability&#45;profile &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/main/artifacts/profile&#45;db)</span> \
                             &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                             &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/main/artifacts/structure&#45;db)</span> \
                             &#45;&#45;engine CDN \
                             &#45;&#45;skip&#45;synonymity
</div>

### Filtering the output 

You can filter the output in various ways, so that you can get straight to the variability positions that you're most interested in. Here are some of the filters that you can set:

* The maximum number of variable positions that can come from a single split (e.g. to look at a max of only two random SCVs from each split)
* The maximum and minimum departure from the reference or consensus position
* The minimum coverage value in all samples (if a position is covered less than that value in a even single sample, it will not be reported)

### Adding additional information

You can also set `--quince-mode`, which reports the variability data across all samples for each position reported (even if that position isn't variable in some samples). For example, if nucleotide position 34 of contig 1 was a SNV in one sample, the output would contain the data for nucleotide position 34 for all of your samples. 

You can also ask the program to report the contig names, split names, and gene-level coverage statistics. 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-variability-profile.md) to update this information.


## Additional Resources


* [All about SNVs, SCVs, and SAAVs](http://merenlab.org/2015/07/20/analyzing-variability/)

* [This program in action in the anvi-3dev tutorial](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#supplying-anvi-3dev-with-sequence-variability)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-variability-profile) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
