---
layout: page
title: anvi-split [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Split an anvi&#x27;o pan or profile database into smaller, self-contained pieces. Provide either a genomes-storage and pan database or a profile and contigs database pair, and you&#x27;ll get back directories of individual projects for each bin  that can be treated as smaller anvi&#x27;o projects.

See **[program help menu](../../../vignette#anvi-split)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[split-bins](../../artifacts/split-bins)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db)</span> <span class="artifact-r">[pan-db](../../artifacts/pan-db)</span> <span class="artifact-r">[collection](../../artifacts/collection)</span></p>

## Usage


This program **creates smaller, self-contained anvi'o projects for each of the <span class="artifact-n">[bin](/software/anvio/help/artifacts/bin)</span>s in your project.** This is useful if you would like to share a subset of an anvi'o project. 

Simply provide either a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> pair or a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/artifacts/genomes-storage-db)</span> and <span class="artifact-n">[pan-db](/software/anvio/help/artifacts/pan-db)</span> pair, as well as a <span class="artifact-n">[collection](/software/anvio/help/artifacts/collection)</span>, and it will create directories for each of your bins that contain their own databases and information. In other words, each of these directories will contain their own anvi'o projects that represent the contigs or genomes stored in that single bin. 

### An example run 

For example, let's say a <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> has a <span class="artifact-n">[collection](/software/anvio/help/artifacts/collection)</span> with three bins, which are (very creatively) called `BIN_1`, `BIN_2`, and `BIN_3`.  

If you ran the following code: 

<div class="codeblock" markdown="1">
anvi&#45;split &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/artifacts/profile&#45;db)</span> \
           &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
           &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/artifacts/collection)</span> \
           &#45;o MY_PATH
</div>

Then in the location `MY_PATH`, you would have three folders: `BIN_1`, `BIN_2`, and `BIN_3`.  Each one contains its own <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> that only contains the contigs from that bin. You can then give a fellow anvi'o user just the `BIN_1` directory and they can get to work. 

Similarly, if you provide a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/artifacts/genomes-storage-db)</span> and <span class="artifact-n">[pan-db](/software/anvio/help/artifacts/pan-db)</span> pair, the directories will contain their own smaller <span class="artifact-n">[genomes-storage-db](/software/anvio/help/artifacts/genomes-storage-db)</span> and <span class="artifact-n">[pan-db](/software/anvio/help/artifacts/pan-db)</span> pairs. 

### Other options 

You are also able to skip generating variability tables or compress the auxiliary data to save space. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-split.md) to update this information.


## Additional Resources


* [Anvi-split in action in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#splitting-the-pangenome)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-split) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
