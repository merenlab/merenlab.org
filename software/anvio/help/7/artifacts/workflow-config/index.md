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

A `JSON`-formated configuration file that describes steps and parameters to be considered by an anvio workflow, which includes <span class="artifact-n">[contigs-workflow](/software/anvio/help/7/artifacts/contigs-workflow)</span>, <span class="artifact-n">[metagenomics-workflow](/software/anvio/help/7/artifacts/metagenomics-workflow)</span>, <span class="artifact-n">[pangenomics-workflow](/software/anvio/help/7/artifacts/pangenomics-workflow)</span>, <span class="artifact-n">[phylogenomics-workflow](/software/anvio/help/7/artifacts/phylogenomics-workflow)</span>, and <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7/artifacts/trnaseq-workflow)</span>.

You can create a default config file for a given workflow using the following command:

```
anvi-run-workflow --workflow ANVIO-WORKFLOW \
                  --get-default-config CONFIG.json
```

Following this, the file `CONFIG.json` will contain all configurable flags and parameters set to their default value for that workflow. From there, you can edit this file to your hearts content. 

### What's in this file? 

The config file contains three types of information:

1. **General parameters**, including the name of the workflow, the version of this config file, and links to the <span class="artifact-n">[fasta-txt](/software/anvio/help/7/artifacts/fasta-txt)</span> or <span class="artifact-n">[samples-txt](/software/anvio/help/7/artifacts/samples-txt)</span> file) 
2. **Rule specific parameters** which allow you to set the parameters on individual anvi'o programs that are run in the workflow. 
3. **Output directory names** which just tell anvi'o what to name all of the intermediate and final outputs (to help keep things organized). 

For example, the default config file for the '<span class="artifact-n">[contigs-workflow](/software/anvio/help/7/artifacts/contigs-workflow)</span>'s no rule specific parameters and looks like this: 

    {
        "workflow_name": "contigs",
        "config_version": 1,
        "fasta_txt": "fasta.txt",
        "output_dirs": {
            "FASTA_DIR":   "01_FASTA_contigs_workflow",
            "CONTIGS_DIR": "02_CONTIGS_contigs_workflow",
            "LOGS_DIR":    "00_LOGS_contigs_workflow"
        }
    }

On the other hand, the default config file for the <span class="artifact-n">[metagenomics-workflow](/software/anvio/help/7/artifacts/metagenomics-workflow)</span> is much longer, because it has sections for each rule specific parameter. For example, its section on parameters for the program <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7/programs/anvi-gen-contigs-database)</span> looks like this:

    "anvi_gen_contigs_database": {
       "--project-name": "{group}",
       "threads": 5,
       "--description": "",
       "--skip-gene-calling": "",
       "--ignore-internal-stop-codons": "",
       "--skip-mindful-splitting": "",
       "--contigs-fasta": "",
       "--split-length": "",
       "--kmer-size": ""
    },

Note that the empty string `""` here means that the default parameter for the program <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7/programs/anvi-gen-contigs-database)</span> will be used. 

For more details on the anvi'o snakemake workflows, please refer to [this tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/).



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/workflow-config.md) to update this information.

