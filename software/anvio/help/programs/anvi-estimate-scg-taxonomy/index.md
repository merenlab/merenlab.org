---
layout: page 
title: anvi-estimate-scg-taxonomy [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/PROGRAM.png" alt="PROGRAM" style="width:100px; border:none" />

Estimates taxonomy at genome and metagenome level. This program is the entry point to estimate taxonomy for a given set of contigs (i.e., all contigs in a contigs database, or contigs described in collections as bins). For this, it uses single-copy core gene sequences and the GTDB database. 

[Back to help main page](../../)

## Provides

<p style="text-align: left" markdown="1"><span style="background:#cbe4d5; padding: 0px 3px 2px 3px; border-radius: 5px;">[genome-taxonomy](../../artifacts/genome-taxonomy)</span> <span style="background:#cbe4d5; padding: 0px 3px 2px 3px; border-radius: 5px;">[genome-taxonomy-txt](../../artifacts/genome-taxonomy-txt)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[profile-db](../../artifacts/profile-db)</span> <span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[contigs-db](../../artifacts/contigs-db)</span> <span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[scgs-taxonomy](../../artifacts/scgs-taxonomy)</span> <span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[collection](../../artifacts/collection)</span> <span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[bin](../../artifacts/bin)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see usage examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs), and feel free to add a Markdown formatted file in this directory named "anvi-estimate-scg-taxonomy.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [Usage examples and warnings](http://merenlab.org/scg-taxonomy)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-estimate-scg-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
