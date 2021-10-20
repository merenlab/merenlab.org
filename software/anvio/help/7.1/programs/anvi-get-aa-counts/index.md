---
layout: page
title: anvi-get-aa-counts [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-get-aa-counts
image:
  featurerelative: ../../../images/header.png
  display: true
---

Fetches the number of times each amino acid occurs from a contigs database in a given bin, set of contigs, or set of genes.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[splits-txt](../../artifacts/splits-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[aa-frequencies-txt](../../artifacts/aa-frequencies-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


Similarly to <span class="artifact-n">[anvi-get-codon-frequencies](/software/anvio/help/7.1/programs/anvi-get-codon-frequencies)</span>, this program counts the number of times each amino acid occurs in a given sequence, whether that's a <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>, <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>, set of contigs (listed in a <span class="artifact-n">[splits-txt](/software/anvio/help/7.1/artifacts/splits-txt)</span>), or a set of genes. The output of this is a <span class="artifact-n">[aa-frequencies-txt](/software/anvio/help/7.1/artifacts/aa-frequencies-txt)</span>. 

There are four possible things you can count the amino acid frequencies in: 
* All of the contigs in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>
* A series of <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>s
* A list of contigs
* A list of genes

Examples for each are below.

### Option 1: all contigs in a contigs-db

To count the amino acids in all of the contigs in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, you can just provide the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> of interest, as so:

<div class="codeblock" markdown="1">
anvi&#45;get&#45;aa&#45;counts &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;o path/to/<span class="artifact&#45;n">[aa&#45;frequencies&#45;txt](/software/anvio/help/7.1/artifacts/aa&#45;frequencies&#45;txt)</span>
</div>

### Option 2: a series of bins in a collection 

To count the amino acid frequencies for a series of <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>s, you'll need to provide three additional parameters: the <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> that you used for binning, the <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> that your bins are contained in, and a text file that describes which bins you are interested in. This text file should have only one bin ID per line. 

So, your run would look something like this: 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;aa&#45;counts &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;o path/to/<span class="artifact&#45;n">[aa&#45;frequencies&#45;txt](/software/anvio/help/7.1/artifacts/aa&#45;frequencies&#45;txt)</span> \
                   &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                   &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> \
                   &#45;B my_favorite_bins.txt
</div>

`my_favorite_bins.txt` would look something like this:

    bin_00001
    bin_00004
    
### Option 3: a list of contigs

Just provide a <span class="artifact-n">[splits-txt](/software/anvio/help/7.1/artifacts/splits-txt)</span> file that lists the contigs you want to look at. 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;aa&#45;counts &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;o path/to/<span class="artifact&#45;n">[aa&#45;frequencies&#45;txt](/software/anvio/help/7.1/artifacts/aa&#45;frequencies&#45;txt)</span> \
                   &#45;&#45;contigs&#45;of&#45;interest <span class="artifact&#45;n">[splits&#45;txt](/software/anvio/help/7.1/artifacts/splits&#45;txt)</span>
</div>

### Option 4: a list of genes 

Just provide a list of gene caller ids, straight into the terminal, like so:

<div class="codeblock" markdown="1">
anvi&#45;get&#45;aa&#45;counts &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;o path/to/<span class="artifact&#45;n">[aa&#45;frequencies&#45;txt](/software/anvio/help/7.1/artifacts/aa&#45;frequencies&#45;txt)</span> \
                   &#45;&#45;gene&#45;caller&#45;ids gene_1,gene_2,gene_3
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-aa-counts.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-aa-counts) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
