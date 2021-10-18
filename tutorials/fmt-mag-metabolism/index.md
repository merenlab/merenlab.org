---
layout: page
title: Metabolic enrichment of high-fitness MAGs from a longitudinal FMT study
modified: 2021-10-22
authors: [iva]
excerpt: "A mini-tutorial on how to run metabolism estimation and enrichment in anvi'o"
categories: [anvio]
comments: true
redirect_from: /fmt-metabolism/
---

This is a **mini-tutorial for the metabolism suite of programs** in anvi'o. We'll be going through how to estimate metabolism and compute enrichment scores for metabolic modules in two different groups of MAGs.

The data we'll be using for this is a real dataset from one of our recent studies, ["Metabolic competency drives microbial colonization and resilience in health and disease”](https://doi.org/10.1101/2021.03.02.433653) by Watson et al. In fact, this post is doing double-duty as a reproducible workflow for one of the analyses in that study. :) The rest of the reproducible workflow can be found [here](https://merenlab.org/data/fmt-gut-colonization/) for anyone who is interested in how we did the other analyses discussed in the paper.

{:.notice}
This tutorial is tailored for anvi'o `v7.1` or later. You can learn the version of your installation by running `anvi-interactive -v` in your terminal.
