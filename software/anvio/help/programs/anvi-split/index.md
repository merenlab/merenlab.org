---
layout: page 
title: anvi-split [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/PROGRAM.png" alt="PROGRAM" style="width:100px; border:none" />

Split an anvi&#39;o pan or profile database into smaller, self-contained pieces. This is usually great when you want to share a subset of an anvi&#39;o project. You give this guy your databases, and a collection id, and it gives you back directories of individual projects for each bin that can be treated as self-contained smaller anvi&#39;o projects. We know you don&#39;t read this far into these help menus, but please remember: you will either need to provide a profile &amp; contigs database pair, or a pan &amp; genomes storage pair. The rest will be taken care of. Magic.

[Back to help main page](../../)

## Provides

<p style="text-align: left" markdown="1"><span style="background:#cbe4d5; padding: 0px 3px 2px 3px; border-radius: 5px;">[split-bins](../../artifacts/split-bins)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[profile-db](../../artifacts/profile-db)</span> <span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[contigs-db](../../artifacts/contigs-db)</span> <span style="background:#dcbfe8; padding: 0px 3px 2px 3px; border-radius: 5px;">[collection](../../artifacts/collection)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see usage examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs), and feel free to add a Markdown formatted file in this directory named "anvi-split.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-split) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
