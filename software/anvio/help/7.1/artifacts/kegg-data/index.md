---
layout: page
title: kegg-data [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/kegg-data
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-setup-kegg-kofams](../../programs/anvi-setup-kegg-kofams)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-run-kegg-kofams](../../programs/anvi-run-kegg-kofams)</span></p>


## Description

A **directory of data** downloaded from the [KEGG database resource](https://www.kegg.jp/) for use in function annotation and metabolism estimation.

It is created by running the program <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7.1/programs/anvi-setup-kegg-kofams)</span>. Not everything from KEGG is included in this directory, only the information relevant to downstream programs. The most critical components of this directory are KOfam HMM profiles and the <span class="artifact-n">[modules-db](/software/anvio/help/7.1/artifacts/modules-db)</span> which contains information on metabolic pathways as described in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html).

Programs that rely on this data directory include <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7.1/programs/anvi-run-kegg-kofams)</span> and <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7.1/programs/anvi-estimate-metabolism)</span>.

## Directory Location
The default location of this data is in the anvi'o folder, at `anvio/anvio/data/misc/KEGG/`. 

You can change this location when you run <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7.1/programs/anvi-setup-kegg-kofams)</span> by providing a different path to the `--kegg-data-dir` parameter:

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;kegg&#45;kofams &#45;&#45;kegg&#45;data&#45;dir /path/to/directory/KEGG
</div>

If you do this, you will need to provide this path to downstream programs that require this data as well.

## Directory Contents

Here is a schematic of how the <span class="artifact-n">[kegg-data](/software/anvio/help/7.1/artifacts/kegg-data)</span> folder will look after setup:

```
KEGG
 |- MODULES.db
 |- ko_list.txt
 |- modules.keg
 |- HMMs
 |   |- Kofam.hmm
 |   |- Kofam.hmm.h3f
 |   |- (....)
 |
 |- modules
 |   |- M00001
 |   |- M00002
 |   |- (....)
 |
 |- orphan_data
     |- 01_ko_fams_with_no_threshold.txt
     |- 02_hmm_profiles_with_ko_fams_with_no_threshold.hmm

```

Typically, users will not have to work directly with any of these files, as downstream programs will interface directly with the <span class="artifact-n">[modules-db](/software/anvio/help/7.1/artifacts/modules-db)</span>. 

However, for the curious:
`ko_list.txt`, `modules.keg`, and all files in the `modules` subfolder are flat text files downloaded from the [KEGG website](https://www.genome.jp/kegg/). The data in these files are processed and organized into the <span class="artifact-n">[modules-db](/software/anvio/help/7.1/artifacts/modules-db)</span> for easier programmatic access. 

The `HMMs` subfolder contains a file of concatentated KOfam profiles (also originally downloaded from [KEGG](https://www.genome.jp/ftp/db/kofam/)), as well as the indexes for this file. Some KOfam profiles do not have a score threshold in the `ko_list.txt` file - these profiles and their corresponding entries from that file live in the `orphan_data` directory. Please note that KOs from the `orphan_data` directory will *not* be annotated in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> when you run <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7.1/programs/anvi-run-kegg-kofams)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/kegg-data.md) to update this information.

