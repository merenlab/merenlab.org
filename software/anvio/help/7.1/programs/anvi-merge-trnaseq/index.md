---
layout: page
title: anvi-merge-trnaseq [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-merge-trnaseq
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program processes one or more anvi&#x27;o tRNA-seq databases produced by `anvi-trnaseq` and outputs anvi&#x27;o contigs and merged profile databases accessible to other tools in the anvi&#x27;o ecosystem. Final tRNA &quot;seed sequences&quot; are determined from a set of samples. Each sample yields a set of tRNA predictions stored in a tRNA-seq database, and these tRNAs may be shared among the samples. tRNA may be 3&#x27; fragments and thereby subsequences of longer tRNAs from other samples which would become seeds. The profile database produced by this program records the coverages of seeds in each sample. This program finalizes predicted nucleotide modification sites using tunable substitution rate parameters..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/semiller10.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Samuel Miller</span><div class="page-author-social-box"><a href="https://semiller10.github.io" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:samuelmiller@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/smiller_science" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/semiller10" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[trnaseq-db](../../artifacts/trnaseq-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[trnaseq-contigs-db](../../artifacts/trnaseq-contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[trnaseq-profile-db](../../artifacts/trnaseq-profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **finds tRNA seed sequences from a set of tRNA-seq samples**.

This program follows <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span> in the <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7.1/artifacts/trnaseq-workflow)</span>. <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span> is run on each tRNA-seq sample, producing sample <span class="artifact-n">[trnaseq-db](/software/anvio/help/7.1/artifacts/trnaseq-db)</span>s. A tRNA-seq database contains predictions of tRNA sequences, structures, and modification sites in the sample. anvi-merge-trnaseq takes as input the tRNA-seq databases from a set of samples. It compares tRNAs predicted from the samples, finding those in common and calculating their sample coverages. The final tRNA sequences predicted from all samples are called **tRNA seeds** and function like contigs in metagenomic experiments. Seeds are stored in a <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span> and sample coverages are stored in a <span class="artifact-n">[trnaseq-profile-db](/software/anvio/help/7.1/artifacts/trnaseq-profile-db)</span>. These databases are **variants** of normal <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>s and <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>s, performing similar functions in the anvi'o ecosystem but containing somewhat different information.

Most of the heavy computational work in the <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7.1/artifacts/trnaseq-workflow)</span> is performed by <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span>. anvi-merge-trnaseq is meant run relatively quickly, allowing its parameters to be tuned to fit the dataset.

The `anvi-merge-trnaseq --help` menu provides detailed explanations of the parameters controlling the multifacted analyses performed by the program.

## Key parameters

### Number of reported seeds

One key parameter is the number of reported tRNA seed sequences (`--max-reported-trna-seeds`). The default value of 10,000 seeds is more appropriate for a complex microbial community than a pure culture of a bacterial isolate, which should yield a number of tRNA seeds equal to the number of expressed tRNAs, say ~30. Sequence artifacts may be reported in addition to the 30 actual tRNAs with a higher value like 10,000. Artifacts are relatively common despite intensive screening by <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span> and anvi-merge-trnaseq due to nontemplated nucleotides and modification-induced mutations introduced into tRNA-seq reads by reverse transcription. In practice, artifacts are easy to distinguish from true tRNA seeds by analyzing seed coverage in <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> and checking seed homology to reference databases, among other measures.

### Modification filters

Other key parameters, `--min-variation` and `--min-third-fourth-nt`, determine the coverage cutoffs that distinguish predicted positions of modified nucleotides from single nucleotide variants. Compared to SNVs, modifications typically produce higher nucleotide variability to three or four different nucleotides. However, modification-induced mutations are often highly skewed to one other nucleotide rather than all three mutant nucleotides. Furthermore, the high coverage of seeds in many tRNA-seq libraries can uncover SNVs with a low-frequency third nucleotide rather than the expected two. Some SNVs that are wrongly called modifications can be easily spotted in <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> and the output of <span class="artifact-n">[anvi-plot-trnaseq](/software/anvio/help/7.1/programs/anvi-plot-trnaseq)</span> due to covariation at two positions in the seed as a result of base pairing. In other words, SNV frequencies are equivalent at the two base paired positions in every sample, where modification artifacts have no effect on nucleotide variability at another position across the molecule.

## Examples

*Merge two samples.*

<div class="codeblock" markdown="1">
anvi&#45;merge&#45;trnaseq trnaseq_database_1 trnaseq_database_2 (...) \
                   &#45;o OUTPUT_DIRECTORY \
                   &#45;n PROJECT_NAME \
</div>

*Merge two samples with and without demethylase treatment, giving priority to the demethylase split in calling the underlying nucleotide at modified positions.*

<div class="codeblock" markdown="1">
anvi&#45;merge&#45;trnaseq untreated_trnaseq_database demethylase_trnaseq_database (...) \
                   &#45;o OUTPUT_DIRECTORY \
                   &#45;n PROJECT_NAME \
                   &#45;&#45;preferred&#45;treatment demethylase
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-merge-trnaseq.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-merge-trnaseq) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
