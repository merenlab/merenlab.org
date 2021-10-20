---
layout: page
title: anvi-oligotype-linkmers [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-oligotype-linkmers
image:
  featurerelative: ../../../images/header.png
  display: true
---

Takes an anvi&#x27;o linkmers report, generates an oligotyping output.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[linkmers-txt](../../artifacts/linkmers-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[oligotypes](../../artifacts/oligotypes) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program converts a <span class="artifact-n">[linkmers-txt](/software/anvio/help/7.1/artifacts/linkmers-txt)</span> artifact into <span class="artifact-n">[oligotypes](/software/anvio/help/7.1/artifacts/oligotypes)</span> data.

A <span class="artifact-n">[linkmers-txt](/software/anvio/help/7.1/artifacts/linkmers-txt)</span> artifact describes each of your short reads that mapped to specific target nucleotide positions in a reference contig. This program counts the total occurance of each combination in those target positions within each of your samples. 

For example, if your <span class="artifact-n">[linkmers-txt](/software/anvio/help/7.1/artifacts/linkmers-txt)</span> focused on two target positions, and you ran the following:

<div class="codeblock" markdown="1">
anvi&#45;oligotype&#45;linkmers &#45;i <span class="artifact&#45;n">[linkmers&#45;txt](/software/anvio/help/7.1/artifacts/linkmers&#45;txt)</span> 
</div>

The output (which by default is called `oligotype-counts-001.txt`) might look like the following:

    key         AG   CA    CG    GA    GG    TA    TG   
    sample_001  0    320   12    2     0     3     579    
    sample_002  0    142   2     0     2     10    353  
    sample_003  3    404   1     1     0     2     610   
    sample_004  0    209   6     0     1     0     240

Note that combinations with zero reads in every sample are not included. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-oligotype-linkmers.md) to update this information.


## Additional Resources


* [An application of the oligotyping workflow in metagenomics](https://merenlab.org/2015/12/09/musings-over-commamox/#an-application-of-oligotyping-in-the-metagenomic-context-oligotyping-amoc)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-oligotype-linkmers) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
