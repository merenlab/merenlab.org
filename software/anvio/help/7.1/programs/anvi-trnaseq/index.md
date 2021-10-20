---
layout: page
title: anvi-trnaseq [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-trnaseq
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program to process reads from a tRNA-seq dataset to generate an anvi&#x27;o tRNA-seq database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/semiller10.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Samuel Miller</span><div class="page-author-social-box"><a href="https://semiller10.github.io" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:samuelmiller@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/smiller_science" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/semiller10" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[trnaseq-fasta](../../artifacts/trnaseq-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[trnaseq-db](../../artifacts/trnaseq-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **analyzes a tRNA-seq library, generating de novo predictions of tRNA sequences, structures, and modification positions**.

A FASTA file of merged paired-end tRNA-seq reads is required as input. This file is produced by the initial steps of the <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7.1/artifacts/trnaseq-workflow)</span>, in which [Illumina-utils](https://github.com/merenlab/illumina-utils), merges paired-end reads and <span class="artifact-n">[anvi-script-reformat-fasta](/software/anvio/help/7.1/programs/anvi-script-reformat-fasta)</span> creates anvi'o-compliant deflines in the FASTA file.

The primary output of anvi-trnaseq is a <span class="artifact-n">[trnaseq-db](/software/anvio/help/7.1/artifacts/trnaseq-db)</span>. Supplemental outputs are also produced -- an analysis summary, a tabular file of unique sequences not identified as tRNA, an a tabular file of 5' and 3' extensions trimmed off mature tRNA.

The `anvi-trnaseq --help` menu provides detailed explanations of the parameters controlling the multifacted analyses performed by the program.

## Examples

*Generate a <span class="artifact-n">[trnaseq-db](/software/anvio/help/7.1/artifacts/trnaseq-db)</span> from a sample using 16 cores.*

<div class="codeblock" markdown="1">
anvi&#45;trnaseq &#45;f <span class="artifact&#45;n">[trnaseq&#45;fasta](/software/anvio/help/7.1/artifacts/trnaseq&#45;fasta)</span> \
             &#45;S SAMPLE_NAME \
             &#45;o OUTPUT_DIRECTORY \
             &#45;T 16
</div>

*Generate a <span class="artifact-n">[trnaseq-db](/software/anvio/help/7.1/artifacts/trnaseq-db)</span> from a sample flagged as being treated with demethylase. The output directory is overwritten if it already exists.*

<div class="codeblock" markdown="1">
anvi&#45;trnaseq &#45;f <span class="artifact&#45;n">[trnaseq&#45;fasta](/software/anvio/help/7.1/artifacts/trnaseq&#45;fasta)</span> \
             &#45;S SAMPLE_NAME \
             &#45;o OUTPUT_DIRECTORY \
             &#45;T 16 \
             &#45;&#45;treatment demethylase \
             &#45;&#45;overwrite&#45;output&#45;destinations
</div>

## Parameterize tRNA feature profiling

Feature profiling parameters can be modified by the user by in an optional `.ini` file. For example, the user may want a more permissive definition of a tRNA (more false positive identifications of sequences as tRNA, fewer false negative failures to identify sequences as tRNA), increasing the number of unpaired nucleotides allowed in the T stem or increasing the number of unconserved canonical nucleotides allowed in the anticodon loop. Numerous structural parameters like these can be altered.

*Write the `.ini` file to `param.ini`.*

<div class="codeblock" markdown="1">
anvi&#45;trnaseq &#45;&#45;default&#45;feature&#45;param&#45;file PARAM.ini
</div>

*Nicely display the `.ini` defaults that can be written to the file in standard output.*

<div class="codeblock" markdown="1">
anvi&#45;trnaseq &#45;&#45;print&#45;default&#45;feature&#45;params
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-trnaseq.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-trnaseq) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
