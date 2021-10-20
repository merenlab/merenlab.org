---
layout: page
title: anvi-script-gen-short-reads [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-gen-short-reads
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate short reads from contigs. Useful to reconstruct mock data sets from already assembled contigs.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[configuration-ini](../../artifacts/configuration-ini) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[short-reads-fasta](../../artifacts/short-reads-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Usage


This program uses already assembled contigs to create a mock list of short reads. You can then use these short reads to reassemble your data in order to test alternative assembly programs or analysis methods as a positive control. 

Basically, this attempts to undo the assembly and produce a data set that could have been directly received from laboratory sequencing. While the computer's mock short reads won't be perfect, they can be used to make sure your analysis pipeline is working from step 1. 

## Example Usage

This program takes an INI file - a form of text file containing various information. For this program, the example provided in the anvi'o test suite looks like this: 

```ini
[general]
short_read_length = 10
error_rate = 0.05
coverage = 100
contig = CTGTGGTTACGCCACCTTGAGAGATATTAGTCGCGTATTGCATCCGTGCCGACAAATTGCCCAACGCATCGTTCCTTCTCCTAAGTAATTTAACATGCGT
```

Note that this file contains both the contig that you want to break down, and various information about the short reads that you want to create. To run this program, just call 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen&#45;short&#45;reads <span class="artifact&#45;n">[configuration&#45;ini](/software/anvio/help/7.1/artifacts/configuration&#45;ini)</span> \
                            &#45;&#45;output&#45;file&#45;path <span class="artifact&#45;n">[short&#45;reads&#45;fasta](/software/anvio/help/7.1/artifacts/short&#45;reads&#45;fasta)</span>
</div>
    
The resulting FASTA file with short reads will cover the `contig` with short reads that are 10 nts long at 100X coverage. There will also be an error-rate of 0.05, to mimic the sequencing errors you would get from sequencing in the wet lab. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen-short-reads.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-short-reads) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
