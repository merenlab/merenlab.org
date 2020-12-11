---
layout: page
title: binding-frequencies-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-interacdome](../../programs/anvi-run-interacdome)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

When <span class="artifact-n">[anvi-run-interacdome](/software/anvio/help/main/programs/anvi-run-interacdome)</span> is ran, it stores binding frequencies directly into the <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> as <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/main/artifacts/misc-data-amino-acids)</span>. Yet <span class="artifact-n">[anvi-run-interacdome](/software/anvio/help/main/programs/anvi-run-interacdome)</span> also outputs tabular data directly accessible by the user--this data is what is meant by <span class="artifact-n">[binding-frequencies-txt](/software/anvio/help/main/artifacts/binding-frequencies-txt)</span>.

Specifically, <span class="artifact-n">[binding-frequencies-txt](/software/anvio/help/main/artifacts/binding-frequencies-txt)</span> refers to 2 files named `INTERACDOME-match_state_contributors.txt` and `INTERACDOME-domain_hits.txt` (the `INTERACDOME` part can be changed with `-O`). One day, the format of these files will be explicitly in this document. Until then, you can find their ouput formats in [this blogpost](https://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database).  Sorry for the inconvenience!


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/binding-frequencies-txt.md) to update this information.

