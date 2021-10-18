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

## A dataset of high- and low-fitness MAGs

You can download the [datapack](FIXME FIGSHARE LINK) for this tutorial by running the following code:

```bash
FIXME DOWNLOAD CODE
```

This dataset includes 40 MAGs of gut microbes, labeled as either "high-fitness" or "low-fitness" according to their [colonization ability](https://merenlab.org/data/fmt-gut-colonization/#defining-colonization-success-and-failure) and prevalence in healthy gut metagenomes (there are 20 MAGs in each group). You can learn the full details of how and why we got them by reading the [study](https://doi.org/10.1101/2021.03.02.433653), but for the purposes of this mini-tutorial, here is what you need to know about these MAGs:

- they were binned from a co-assembly of **longitudinally-sampled gut metagenomes** taken from a **healthy adult who donated stool for fecal microbiota transplantation (FMT)**
- the **high-fitness MAGs** represent microbial **populations that were able to colonize all FMT recipients** who received stool from this donor. They were detected (with sufficient abundance) in recipient gut metagenomes at least 7 days (and up to 1 year) post-FMT
- the **low-fitness MAGs** represent **populations that were NOT able to colonize** FMT recipients who received stool from this donor. They were not detected in recipient gut metagenomes post-FMT
- the **high-fitness MAGs** were the 20 MAGs with the **highest prevalence in gut metagenomes from healthy Canadian adults**, while the **low-fitness MAGs were less prevalent**

The "high-fitness" MAGs were labeled that way because we hypothesized that something about these populations increased their fitness such that they were able to survive the stress of being transplanted into new gut environments and become long-term colonizers (in comparison to the "low-fitness" populations that were unable to survive for long in the recipients). In our study, we sought to learn what distinguishes these two groups from each other - what enables one group to survive while the other does not? What do the "high-fitness" populations have that the "low-fitness" ones don't (or vice versa)?

One way to answer this question is to look at the metabolic potential, or genomically-encoded metabolic capabilities, of these MAGs.
