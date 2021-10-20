---
layout: page
title: interacdome-data [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/interacdome-data
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DATA.png" alt="DATA" style="width:100px; border:none" />

A DATA-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-setup-interacdome](../../programs/anvi-setup-interacdome)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-run-interacdome](../../programs/anvi-run-interacdome)</span></p>


## Description

This artifact stores the data downloaded by <span class="artifact-n">[anvi-setup-interacdome](/software/anvio/help/7.1/programs/anvi-setup-interacdome)</span> and is required to run <span class="artifact-n">[anvi-run-interacdome](/software/anvio/help/7.1/programs/anvi-run-interacdome)</span>.

As described in the [InteracDome blogpost](https://merenlab.org/2020/07/22/interacdome/#anvi-setup-interacdome), this data includes [the tab-separated files](https://interacdome.princeton.edu/#tab-6136-4) from [InteracDome](https://interacdome.princeton.edu/) and the Pfam 31.0 HMMs for the subset of Pfams found in the InteracDome dataset.

By default, this data is stored in `anvio/data/misc/InteracDome`, but a custom path can be set when the user runs <span class="artifact-n">[anvi-setup-interacdome](/software/anvio/help/7.1/programs/anvi-setup-interacdome)</span> if desired.



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/interacdome-data.md) to update this information.

