---
layout: page
title: anvi-run-kegg-kofams [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-run-kegg-kofams
image:
  featurerelative: ../../../images/header.png
  display: true
---

Run KOfam HMMs on an anvi&#x27;o contigs database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[kegg-data](../../artifacts/kegg-data) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[kegg-functions](../../artifacts/kegg-functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[functions](../../artifacts/functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


Essentially, this program uses the KEGG database to annotate functions and metabolic pathways in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. More specifically, <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7.1/programs/anvi-run-kegg-kofams)</span> annotates a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> with HMM hits from KOfam, a database of KEGG Orthologs (KOs). You must set up these HMMs on your computer using <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7.1/programs/anvi-setup-kegg-kofams)</span> before you can use this program.

Running this program is a pre-requisite for metabolism estimation with <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7.1/programs/anvi-estimate-metabolism)</span>. Note that if you are planning to run metabolism estimation, it must be run with the same <span class="artifact-n">[kegg-data](/software/anvio/help/7.1/artifacts/kegg-data)</span> that is used in this program to annotate KOfam hits.

### How does it work?
**1) Run an HMM search against KOfam**
Briefly, what this program does is extract all the gene calls from the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and checks each one for hits to the KOfam HMM profiles in your <span class="artifact-n">[kegg-data](/software/anvio/help/7.1/artifacts/kegg-data)</span>. This can be time-consuming given that the number of HMM profiles is quite large, even more so if the number of genes in the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> is also large. Multi-threading is a good idea if you have the computational capability to do so.

**2) Eliminate weak hits based on bitscore**
Many HMM hits will be found, most of them weak. The weak hits will by default be eliminated according to the bitscore thresholds provided by KEGG; that is, hits with bitscores below the threshold for a given KO profile will be discarded, and those with bitscores above the threshold will be annotated in the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. It is perfectly normal to notice that the number of raw hits found is many, many times larger than the number of annotated KO hits in your database.

**3) Add back valid hits that were missed**
There is one issue with this practice of removing _all_ KOfam hits below the KEGG bitscore threshold for a given profile. We (and others) have noticed that the KEGG thresholds can sometimes be too stringent, eliminating hits that are actually valid annotations. To solve this problem, we
have implemented the following heuristic for relaxing the bitscore thresholds and annotating genes that would otherwise go without a valid KO annotation:

For every gene without a KOfam annotation, we examine all the hits with an e-value below `X` and a bitscore above `Y` percent of the threshold. If those hits are all to a unique KOfam profile, then we annotate the gene call with that KO.

`X` and `Y` are parameters that can be modified (see below), but by default the e-value threshold (`X`) is 1e-05 and the bitscore fraction (`Y`) is 0.5.

Please note that this strategy is just a heuristic. We have tried to pick default parameters that seemed reasonable but by no means have we comprehensively tested and optimized them. This is why X and Y are mutable so that you can explore different values and see how they work for your data. It is always a good idea to double-check your annotations to make sure they are reasonable and as stringent as you'd like them to be. In addition, if you do not feel comfortable using this heuristic at all, you can always turn this behavior off and rely solely on KEGG's bitscore thresholds. :)

**3) Put annotations in the database**
In the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> functions table, annotated KO hits (<span class="artifact-n">[kegg-functions](/software/anvio/help/7.1/artifacts/kegg-functions)</span>) will have the source `KOfam`.

### Standard usage

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db
</div>

### Use a specific non-default KEGG data directory
If you have previously setup your KEGG data directory using `--kegg-data-dir` (see <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7.1/programs/anvi-setup-kegg-kofams)</span>), or have moved the KEGG data directory that you wish to use to a non-default location (maybe you like keeping the older versions around when you update, we don't know how you roll), then you may need to specify where to find the KEGG data so that this program can use the right one. In that case, this is how you do it:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db \
                     &#45;&#45;kegg&#45;data&#45;dir /path/to/directory/KEGG
</div>

### Run with multiple threads

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db &#45;T 4
</div>

### Use a different HMMER program
By default, <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7.1/programs/anvi-run-kegg-kofams)</span> uses `hmmsearch` to find KO hits. If for some reason you would rather use a different program (`hmmscan` is also currently supported), you can do so.

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db \
                     &#45;&#45;hmmer&#45;program hmmscan
</div>

### Keep all HMM hits
Usually, this program parses out weak HMM hits and keeps only those that are above the score threshold for a given KO. If you would like to turn off this behavior and keep all hits (there will be _a lot_ of weak ones), you can follow the example below:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db \
                     &#45;&#45;keep&#45;all&#45;hits
</div>

### Modifying the bitscore relaxation heuristic
As described above, this program does its best to avoid missing valid annotations by relaxing the bitscore threshold for genes without any annotations. For such a gene, hits with e-value <= X and bitscore > (Y * KEGG threshold) that are all hits to the same KOfam profile are used to annotate the gene with that KO.

**Skip this step entirely**
If you don't want any previously-eliminated hits to be used for annotation, you can skip this heuristic by using the flag `--skip-bitscore-heuristic`. Then, _only_ hits with bitscores above the KEGG-provided threshold for a given KO will be used for annotation.

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db \
                     &#45;&#45;skip&#45;bitscore&#45;heuristic
</div>

**Modify the heuristic parameters**
If our default values are too stringent or not stringent enough for your tastes, you can change them! The e-value threshold (X, default: 1e-05) can be set using `-E` or `--heuristic-e-value` and the bitscore fraction (Y, default: 0.50) can be set using `-H` or `--heuristic-bitscore-fraction`. Like so:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db \
                     &#45;E 1e&#45;15 \
                     &#45;H 0.90
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-kegg-kofams.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-kegg-kofams) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
