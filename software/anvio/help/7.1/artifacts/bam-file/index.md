---
layout: page
title: bam-file [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/bam-file
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/BAM.png" alt="BAM" style="width:100px; border:none" />

A BAM-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-init-bam](../../programs/anvi-init-bam)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-get-short-reads-from-bam](../../programs/anvi-get-short-reads-from-bam)</span> <span class="artifact-r">[anvi-get-short-reads-mapping-to-a-gene](../../programs/anvi-get-short-reads-mapping-to-a-gene)</span> <span class="artifact-r">[anvi-get-tlen-dist-from-bam](../../programs/anvi-get-tlen-dist-from-bam)</span> <span class="artifact-r">[anvi-profile](../../programs/anvi-profile)</span> <span class="artifact-r">[anvi-profile-blitz](../../programs/anvi-profile-blitz)</span> <span class="artifact-r">[anvi-report-linkmers](../../programs/anvi-report-linkmers)</span> <span class="artifact-r">[anvi-script-get-coverage-from-bam](../../programs/anvi-script-get-coverage-from-bam)</span></p>


## Description

A BAM file contains **already aligned sequence data.** However, it is written in binary to save space (so it will look like jibberish if you open it). 

BAM files (and their text file cousin SAM files) are often used in 'omics analysis and are described in more detail in [this file](https://samtools.github.io/hts-specs/SAMv1.pdf), written by the developers of samtools. 

If your BAM file is not indexed, it is actually a <span class="artifact-n">[raw-bam-file](/software/anvio/help/7.1/artifacts/raw-bam-file)</span> and you can run <span class="artifact-n">[anvi-init-bam](/software/anvio/help/7.1/programs/anvi-init-bam)</span> to turn it into a BAM file. You can tell if your BAM file is indexed if in the same folder as your `XXXX.bam` file, there is another file with the same name called `XXXX.bam.bai`.

As of now, no anvi'o programs will output results in BAM format, so you'll primary use BAM files to import sequence data into anvi'o. For example, in <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span> (which generates a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>), your BAM file is expected to contain the aligned short reads from your samples. 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/bam-file.md) to update this information.

