---
layout: page
title: anvi-script-augustus-output-to-external-gene-calls [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-augustus-output-to-external-gene-calls
image:
  featurerelative: ../../../images/header.png
  display: true
---

Takes in gene calls by AUGUSTUS v3.3.3, generates an anvi&#x27;o external gene calls file. It may work well with other versions of AUGUSTUS, too. It is just no one has tested the script with different versions of the program.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[augustus-gene-calls](../../artifacts/augustus-gene-calls) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[external-gene-calls](../../artifacts/external-gene-calls) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program converts a gene call file from [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) (as an <span class="artifact-n">[augustus-gene-calls](/software/anvio/help/7.1/artifacts/augustus-gene-calls)</span> artifact) to an anvi'o <span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span> artifact. 

This essentially just reformats the data in the <span class="artifact-n">[augustus-gene-calls](/software/anvio/help/7.1/artifacts/augustus-gene-calls)</span> artifact (for example, removing the UTR information) so that it can be read by other anvi'o programs. 

A run of this program will look something like this:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;augustus&#45;output&#45;to&#45;external&#45;gene&#45;calls &#45;i <span class="artifact&#45;n">[augustus&#45;gene&#45;calls](/software/anvio/help/7.1/artifacts/augustus&#45;gene&#45;calls)</span>
                                                   &#45;o <span class="artifact&#45;n">[external&#45;gene&#45;calls](/software/anvio/help/7.1/artifacts/external&#45;gene&#45;calls)</span>
</div>

Here is an example of what the resulting <span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span> file will look like (from the gff file used as an example on the <span class="artifact-n">[augustus-gene-calls](/software/anvio/help/7.1/artifacts/augustus-gene-calls)</span> page):  

    gene_callers_id    contig       start    stop    direction    partial    call_type    source      version    aa_sequence
    0                  unnamed-1    56       1252    f            0          1            AUGUSTUS    v3.3.3     MSEGNAAGEPSTPGGPRPLLTGARGLIGRRPAPPLTPGRLPSIRSRDLTLGGVKKKTFTPNIISRKIKEEPKEEVTVKKEKRERDRDRQREGHGRGRGRPEVIQSHSIFEQGPAEMMKKKGNWDKTVDVSDMGPSHIINIKKEKRETDEETKQILRMLEKDDFLDDPGLRNDTRNMPVQLPLAHSGWLFKEENDEPDVKPWLAGPKEEDMEVDIPAVKVKEEPRDEEEEAKMKAPPKAARKTPGLPKDVSVAELLRELSLTKEEELLFLQLPDTLPGQPPTQDIKPIKTEVQGEDGQVVLIKQEKDREAKLAENACTLADLTEGQVGKLLIRKSGRVQLLLGKVTLDVTMGTACSFLQELVSVGLGDSRTGEMTVLGHVKHKLVCSPDFESLLDHKHR



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-augustus-output-to-external-gene-calls.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-augustus-output-to-external-gene-calls) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
