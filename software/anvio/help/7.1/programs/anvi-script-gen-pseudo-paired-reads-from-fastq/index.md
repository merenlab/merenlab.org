---
layout: page
title: anvi-script-gen-pseudo-paired-reads-from-fastq [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-gen-pseudo-paired-reads-from-fastq
image:
  featurerelative: ../../../images/header.png
  display: true
---

A script that takes a FASTQ file that is not paired-end (i.e., R1 alone) and converts it into two FASTQ files that are paired-end (i.e., R1 and R2). This is a quick-and-dirty workaround that halves each read from the original FASTQ and puts one half in the FASTQ file for R1 and puts the reverse-complement of the second half in the FASTQ file for R2. If you&#x27;ve ended up here, things have clearly not gone very well for you, and Evan, who battled similar battles and ended up implementing this solution wholeheartedly sympathizes.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[short-reads-fasta](../../artifacts/short-reads-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[short-reads-fasta](../../artifacts/short-reads-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Usage


This program takes in a <span class="artifact-n">[short-reads-fasta](/software/anvio/help/7.1/artifacts/short-reads-fasta)</span> file and tries to recreate what paired reads for the data in that fasta file might look like. 

An arbitrarily chosen half of the reads will be put into the R1 output, while the other half will be reverse complemented and put into the R2 output. 

For example, if you ran 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen&#45;pseudo&#45;paired&#45;reads&#45;from&#45;fastq &#45;f <span class="artifact&#45;n">[short&#45;reads&#45;fasta](/software/anvio/help/7.1/artifacts/short&#45;reads&#45;fasta)</span> \
                                               &#45;O MY_READS 
</div>

Then you would end up with two files: 

- `MY_READS_1.fastq` which contains half of the reads straight out of your input file
- `MY_READS_2.fastq` which contains the reverse complement of the other half of the reads. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen-pseudo-paired-reads-from-fastq.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-pseudo-paired-reads-from-fastq) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
