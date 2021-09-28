---
layout: post
authors: [iva]
title: "Targeted binning of a novel nitrogen-fixing population from the Arctic Ocean"
excerpt: "Using metabolism estimation to go fishing in some ocean metagenomes."
modified: 2021-10-20
tags: [metabolism]
categories: [miscellaneous, anvio]
comments: true
redirect_from:
  - /targeted-binning-nif-mag/
---

Public sequencing datasets are beautiful things. They open a world of endless possibilities to scientists everywhere simply by being downloadable. It doesn't matter what they were originally sequenced for, because curious people can always find something to explore in them. Stool metagenomes from microbiome studies like [this one](https://doi.org/10.1016/j.cub.2015.04.055) might eventually be used [to find human gut-associated plasmids](https://www.biorxiv.org/content/10.1101/2020.11.01.361691v1.full). Genomes and MAGs from sequencing projects all over the globe come together in the phylogenies and pangenomes of countless papers. Today's [time-series](https://pubmed.ncbi.nlm.nih.gov/22936250/) can turn into tomorrow's [tutorial](https://merenlab.org/tutorials/infant-gut/). That's the magic of open science.

At the same time, the curse of open science is that sequencing projects designed and carried out by other people are not often obviously or directly usable in your own projects. Usually, you don't find what you are looking for even if you have the right tools to look for it. But sometimes, you get lucky. This post is about one of those times.

If you are interested in learning about how to leverage anvi'o's metabolism estimation capabilities to go fishing through your data, or if you just get a kick out of success stories in the realm of open science and data reusability, keep on reading.

{:.notice}
This post also doubles as a reproducible workflow. Feel free to download the associated datapack at [INSERT DATAPACK LINK]() and follow along with the commands. (Or, go your own way and explore the data yourself!) The commands below were written for anvi'o `v7`. If you have a newer version of anvi'o and find a command that is not working - sorry. We don't always retroactively update these posts as anvi'o evolves. Please send us a message and we'll see if we can fix it.  

(TODO: check v7 compatibility of commands)

## Nitrogen fixation is a big thing
... in the temperate oceans, at least.
