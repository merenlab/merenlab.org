---
layout: page
title: workflow-config [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/JSON.png" alt="JSON" style="width:100px; border:none" />

A JSON-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-migrate](../../programs/anvi-migrate)</span> <span class="artifact-r">[anvi-run-workflow](../../programs/anvi-run-workflow)</span></p>


## Description

A `JSON`-formated configuration file that describes steps and parameters to be considered by an anvio workflow, which includes <span class="artifact-n">[contigs-workflow](/software/anvio/help/main/artifacts/contigs-workflow)</span>, <span class="artifact-n">[metagenomics-workflow](/software/anvio/help/main/artifacts/metagenomics-workflow)</span>, <span class="artifact-n">[pangenomics-workflow](/software/anvio/help/main/artifacts/pangenomics-workflow)</span>, <span class="artifact-n">[phylogenomics-workflow](/software/anvio/help/main/artifacts/phylogenomics-workflow)</span>, and <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/main/artifacts/trnaseq-workflow)</span>.

You can create a default config file for a gein workflow using the following command:

```
anvi-run-workflow --workflow ANVIO-WORKFLOW \
                  --get-default-config CONFIG.json
```

For details of anvi'o snakemake workflows, please refer to [this tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/workflow-config.md) to update this information.

