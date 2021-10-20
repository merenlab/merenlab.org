---
layout: page
title: anvi-run-hmms [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-run-hmms
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program deals with populating tables that store HMM hits in an anvi&#x27;o contigs database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/mschecht.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Matthew Schechter</span><div class="page-author-social-box"><a href="mailto:mschechter@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/mschecht_bio" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/mschecht" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-source](../../artifacts/hmm-source) <img src="../../images/icons/HMM.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[hmm-hits](../../artifacts/hmm-hits) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


Stores <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span> for a given <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span> in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. In short, this is the program that will do a search for HMMs against a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and store that information into the contigs-db's <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span>.

This is one of the programs that users commonly run on newly generated <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, along with <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/7.1/programs/anvi-scan-trnas)</span>, <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/7.1/programs/anvi-run-ncbi-cogs)</span>, <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-run-scg-taxonomy)</span>, and so on.

### HMMs in the context of anvi'o

In a nutshell, [hidden Markov models](https://en.wikipedia.org/wiki/Hidden_Markov_model) are statistical models typically generated from known genes which enable 'searching' for similar genes in other sequence contexts.

The default anvi'o distribution includes numerous [curated HMM profiles](https://github.com/merenlab/anvio/tree/master/anvio/data/hmm) for single-copy core genes and ribosomal RNAs, and anvi'o can work with custom HMM profiles provided by the user. In anvi'o lingo, each of these HMM profiles, whether they are built-in or user defined, is called an <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>.

### Default Usage

To run this program with all default settings (against all default anvi'o <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>), you only need to provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span>
</div>

Multithreading will dramatically improve the performance of `anvi-run-hmms`. If you have multiple CPUs or cores, you may parallelize your search:


<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
              &#45;&#45;num&#45;threads 6
</div>


You can also run this program on a specific built-in <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
              &#45;I Bacteria_71
</div>

### User-defined HMMs

Running `anvi-run-hmms` with a custom model is easy. All you need to do is to create a directory with necessary files:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
              &#45;H MY_HMM_PROFILE
</div>

See the relevant section in the artifact <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span> for details.


### Changing the HMMER program

By default, `anvi-run-hmms` will use [HMMER](http://hmmer.org/)'s `hmmscan` for amino acid HMM profiles, but you can use `hmmsearch` if you are searching a very large number of models against a relatively smaller number of sequences:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
              &#45;&#45;hmmer&#45;program hmmsearch
</div>

{:.notice}
This flag has no effect when your HMM profile source is for nucleotide sequences (like any of the Ribosomal RNA sources). In those cases anvi'o will use `nhmmscan` exclusively.

### Saving the HMMER output

If you want to see the output from the HMMER program (eg, `hmmscan`) used to annotate your data, you can request that it be saved in a directory of your choosing. Please note that this only works when you are running on a single HMM source, as in the example below:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
              &#45;I Bacteria_71 \
              &#45;&#45;hmmer&#45;output&#45;dir OUTPUT_DIR
</div>

If you do this, file(s) with the prefix `hmm` will appear in that directory, with the file extension indicating the format of the output file. For example, the table output format would be called `hmm.table`.

{:.warning}
These resulting files are not _exactly_ the raw output of HMMER because anvi'o does quite a bit of pre-processing on the raw input and output file(s) while jumping through some hoops to make the HMM searches multi-threaded. If this is causing you a lot of headache, please let us know.

#### Requesting domain table output

{:.notice}
Please also see <span class="artifact-n">[anvi-script-filter-hmm-hits-table](/software/anvio/help/7.1/programs/anvi-script-filter-hmm-hits-table)</span>

No matter what, anvi'o will use the regular table output to annotate your contigs database. However, if you are using the `--hmmer-output-dir` to store the HMMER output, you can also request a domain table output using the flag `--domain-hits-table`.

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
              &#45;I Bacteria_71 \
              &#45;&#45;hmmer&#45;output&#45;dir OUTPUT_DIR \
              &#45;&#45;domain&#45;hits&#45;table
</div>

In this case anvi'o will run [HMMER](http://hmmer.org) using the `--domtblout` flag to generate this output file.

{:.notice}
This flag will only work with HMM profiles made for amino acid sequences. Profiles for nucleotide sequences require the use of the program `nhmmscan`, which does not have an option to store domain output.

Please note that this output **won't be used to filter hits to be added to the contigs database**. But it will give you the necessary output file to investigate the coverage of HMM hits. But you can use the program <span class="artifact-n">[anvi-script-filter-hmm-hits-table](/software/anvio/help/7.1/programs/anvi-script-filter-hmm-hits-table)</span> with this file to remove weak hits from your HMM hits table later.


### Other things anvi-run-hmms can do

* Add the tag `--also-scan-trnas` to basically run <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/7.1/programs/anvi-scan-trnas)</span> for you at the same time. It's very convenient. (But it only works if you are not using the `-I` or `-H` flags at the same time because reasons.)



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-hmms.md) to update this information.


## Additional Resources


* [Another description as part of the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-hmms) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
