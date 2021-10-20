---
layout: page
title: anvi-script-filter-fasta-by-blast [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-filter-fasta-by-blast
image:
  featurerelative: ../../../images/header.png
  display: true
---

Filter FASTA file according to BLAST table (remove sequences with bad BLAST alignment).

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/no-avatar.png" /></div><div class="page-person-info-box"><span class="page-author-name">Daniel Blankenberg</span><div class="page-author-social-box"><a href="mailto:blanked2@ccf.org" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/blankenberg" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-fasta](../../artifacts/contigs-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[blast-table](../../artifacts/blast-table) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-fasta](../../artifacts/contigs-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Usage


This program takes a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span> and <span class="artifact-n">[blast-table](/software/anvio/help/7.1/artifacts/blast-table)</span> and removes sequences without BLAST hits of a certain level of confidence. 

For example, you could use this program to filter out sequences that do not have high-confidence taxonomy assignments before running a phylogenomic analysis. 

To run this program, you'll need to provide the <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span> that you're planning to filter, the <span class="artifact-n">[blast-table](/software/anvio/help/7.1/artifacts/blast-table)</span>, a list of the column headers in your <span class="artifact-n">[blast-table](/software/anvio/help/7.1/artifacts/blast-table)</span> (as given to BLAST by `-outfmt`), and a `proper_pident` threshold at which to remove the sequences. This threshold will remove sequences less than the given percent of the query amino acids that were identical to the corresponding matched amino acids. Note that this diffres from the `pident` blast parameter because it doesn't include unaligned regions. 

For example, if you ran 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;filter&#45;fasta&#45;by&#45;blast &#45;f <span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/7.1/artifacts/contigs&#45;fasta)</span> \
                                  &#45;o path/to/<span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/7.1/artifacts/contigs&#45;fasta)</span> \
                                  &#45;b <span class="artifact&#45;n">[blast&#45;table](/software/anvio/help/7.1/artifacts/blast&#45;table)</span> \
                                  &#45;s qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
                                  &#45;t 30
</div>
        
Then the output file would be a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span> that contains only the sequences in your input file that have a hit in your blast table with more than 30 percent of the amino acids aligned. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-filter-fasta-by-blast.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-filter-fasta-by-blast) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
