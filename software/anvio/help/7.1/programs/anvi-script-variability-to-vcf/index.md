---
layout: page
title: anvi-script-variability-to-vcf [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-script-variability-to-vcf
image:
  featurerelative: ../../../images/header.png
  display: true
---

A script to convert SNV output obtained from anvi-gen-variability-profile to the standard VCF format.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/no-avatar.png" /></div><div class="page-person-info-box"><span class="page-author-name">Srinidhi Varadharajan</span><div class="page-author-social-box"><a href="mailto:srinidhi.varadharajan@bcm.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/srinidhi202" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[vcf](../../artifacts/vcf) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This script **converts a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/main/artifacts/variability-profile-txt)</span> into <span class="artifact-n">[vcf](/software/anvio/help/main/artifacts/vcf)</span> (Variant Call Format).** 

It is very easy to run: just provide the input and output paths as so:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;variability&#45;to&#45;vcf &#45;i <span class="artifact&#45;n">[variability&#45;profile&#45;txt](/software/anvio/help/main/artifacts/variability&#45;profile&#45;txt)</span> \ 
                               &#45;o <span class="artifact&#45;n">[vcf](/software/anvio/help/main/artifacts/vcf)</span> 
</div>

Note that to run this, you'll need to have run <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/main/programs/anvi-gen-variability-profile)</span> with the default nucleotide engine. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-variability-to-vcf.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-variability-to-vcf) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
