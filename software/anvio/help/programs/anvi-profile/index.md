---
layout: page 
title: anvi-profile [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Creates a single anvi&#39;o profile database. The default input to this program is a BAM file.                   When it is run on a BAM file, depending on the user parameters, the program quantifies                   coverage per nucleotide position (and averages them out per contig), calculates                   single-nucleotide, single-codon, and single-amino acid variants, and stores these data                   into appropriate tables. Anvi&#39;o single profiles can be merged by the program `anvi-merge`.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


[Back to help main page](../../)

## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[single-profile-db](../../artifacts/single-profile-db)</span> <span class="artifact-p">[misc-data-item-orders](../../artifacts/misc-data-item-orders)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[bam-file](../../artifacts/bam-file)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-profile.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [The usage of the profiler in metagenomics workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-profile) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
