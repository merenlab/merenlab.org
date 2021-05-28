---
layout: page
title: modules-db [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-setup-kegg-kofams](../../programs/anvi-setup-kegg-kofams)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-migrate](../../programs/anvi-migrate)</span></p>


## Description

A database containing information from the [KEGG MODULES database](https://www.genome.jp/kegg/module.html) for use in metabolic reconstruction and functional annotation of KEGG Orthologs (KOs).

Part of the <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span> directory. You can get this database on your computer by running <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7/programs/anvi-setup-kegg-kofams)</span>. Programs that rely on this database downstream include <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span> and <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/modules-db.md) to update this information.

