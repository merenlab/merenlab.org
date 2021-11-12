---
layout: page
title: anvi-export-gene-calls [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-export-gene-calls
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export gene calls from an anvi&#x27;o contigs database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[gene-calls-txt](../../artifacts/gene-calls-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


The purpose of this program is to exports your gene calls in a given <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> and a gene caller, in the form of a <span class="artifact-n">[gene-calls-txt](/software/anvio/help/main/artifacts/gene-calls-txt)</span>. 

To see the gene callers available in your contigs database, you can use <span class="artifact-n">[anvi-db-info](/software/anvio/help/main/programs/anvi-db-info)</span> or use this program with the following flag: 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;gene&#45;calls &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                       &#45;&#45;list&#45;gene&#45;callers
</div>

Running this will export all of your gene calls identified by the gene caller [prodigal](https://github.com/hyattpd/Prodigal) (assuming it is in your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>:

<div class="codeblock" markdown="1">
anvi&#45;export&#45;gene&#45;calls &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                       &#45;&#45;gene&#45;caller Prodigal \
                       &#45;o <span class="artifact&#45;n">[gene&#45;calls&#45;txt](/software/anvio/help/main/artifacts/gene&#45;calls&#45;txt)</span>
</div>

{:.notice}
You can export genes from more gene callers by providing a comma-separated list of gene caller names.

If you don't want to display the amino acid sequences of each gene (they can crowd the file very quickly if you don't want to see them), you can add the following flag:

<div class="codeblock" markdown="1">
anvi&#45;export&#45;gene&#45;calls &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                       &#45;&#45;gene&#45;caller Prodigal \
                       &#45;&#45;skip&#45;sequence&#45;reporting \
                       &#45;o <span class="artifact&#45;n">[gene&#45;calls&#45;txt](/software/anvio/help/main/artifacts/gene&#45;calls&#45;txt)</span>
</div>

## Advanced uses

This program can take a lot of time and memory when working with very large <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> files (such as those that are more than 10 Gb in file size or more than 10 million contigs).

In that case you can export your gene calls the following way within minutes and a small memory space.

First open your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>:

<div class="codeblock" markdown="1">
sqlite3 <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span>
</div>

Then run these lines,

<div class="codeblock" markdown="1">
.mode csv 
.headers on 
.out <span class="artifact&#45;n">[gene&#45;calls&#45;txt](/software/anvio/help/main/artifacts/gene&#45;calls&#45;txt)</span>
select gene_callers_id, contig, start, stop, direction, partial from genes_in_contigs;
</div>

You can also continue with these lines to get the amino acid sequences for them:

<div class="codeblock" markdown="1">
.mode csv 
.headers on 
.out AMINO&#45;ACID&#45;SEQUENCES.txt
select &#42; from genes_in_contigs;
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-gene-calls.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-gene-calls) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
