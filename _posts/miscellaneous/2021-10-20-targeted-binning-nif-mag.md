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
This post also doubles as a reproducible workflow. Feel free to download the associated datapack at [FIXME INSERT DATAPACK LINK]() and follow along with the commands. (Or, go your own way and explore the data yourself!) The commands below were written for anvi'o `v7`. If you have a newer version of anvi'o and find a command that is not working - sorry. We don't always retroactively update these posts as anvi'o evolves. Please send us a message and we'll see if we can fix it.  

(TODO: check v7 compatibility of commands)

## A focus on nitrogen fixation

About one year ago, I was looking for nitrogen fixation genes in ocean metagenomes. To tell you the truth, it wasn't so much that I cared a lot about nitrogen fixation back then (I had only just started getting into the whole metabolism thing), and more that I really just wanted to test out [this new program that I wrote](https://merenlab.org/software/anvio/help/7/programs/anvi-estimate-metabolism/) in as many different metagenome types as I could. So when my much-more-knowledgable colleagues suggested that I use this tool to search for nitrogen-fixing populations in a non-conventional environment - polar ocean metagenomes - I said, "Why not?"

For anyone who is as clueless as I was about nitrogen fixation, here is a very light summary of the background. Nitrogen fixation is the process of converting gaseous nitrogen (N<sub>2</sub>) into the more biologically-usable form ammonia (NH<sub>3</sub>). The resulting ammonia can then be further converted into other bioavailable compounds like nitrate and nitrite. Since nitrogen is an essential component of many biological molecules (amino acids, anyone?), nitrogen fixation is a fairly important process that supports life, in general. Only some ocean microbes, called marine diazotrophs, have the genes that allow them to fix nitrogen, which are encoded in the _nif_ operon. ([Sohm et al 2011](https://www.nature.com/articles/nrmicro2594))

Why are polar oceans a non-conventional place to find nitrogen fixers? Well, the majority of nitrogen fixing microbes have been found in non-polar oceans, perhaps because the diazotrophic cyanobacteria that are typically studied tend to be found in warmer waters ([Stal 2009](https://doi.org/10.1111/j.1758-2229.2009.00016.x)), or perhaps because the polar oceans are just not as well-studied as the other oceans.

To get back to the story, in early 2020 I started looking for nitrogen fixation genes in some private polar ocean metagenomes, graciously shared by our collaborators at FIXME (who will hopefully share their findings in an exciting study in the future). And I found nothing. Nada. Not a single nitrogen fixation gene. I promise I tried very hard. But you cannot find things if they are not there.

However, Meren eventually brought to my attention a recently published dataset of Artic and Antarctic ocean metagenomes.
