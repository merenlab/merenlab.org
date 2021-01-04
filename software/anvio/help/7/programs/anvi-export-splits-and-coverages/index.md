---
layout: page
title: anvi-export-splits-and-coverages [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export split or contig sequences and coverages across samples stored in an anvi&#x27;o profile database. This program is especially useful if you would like to &#x27;bin&#x27; your splits or contigs outside of anvi&#x27;o and import the binning results into anvi&#x27;o using `anvi-import-collection` program.

See **[program help menu](../../../../vignette#anvi-export-splits-and-coverages)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-fasta](../../artifacts/contigs-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[coverages-txt](../../artifacts/coverages-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program **gives you the coverage information in your <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> as external files**. Basically, if you want to take that information in your <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> out of anvio, this is for you. 

Once you input your <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> and the <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> you used to generate it, it will create a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7/artifacts/contigs-fasta)</span> that lists your contigs for you, as well as a <span class="artifact-n">[coverages-txt](/software/anvio/help/7/artifacts/coverages-txt)</span>, which describes your coverage information. 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;splits&#45;and&#45;coverages &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span>
</div>

If your coverages are skewed by outlier positions, consider using Q2Q3-coverages instead.

<div class="codeblock" markdown="1">
anvi&#45;export&#45;splits&#45;and&#45;coverages &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                                 &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                                 &#45;&#45;use&#45;Q2Q3&#45;coverages
</div>

### Contigs or splits?

*Wondering what the difference is? Check out [our vocab page](http://merenlab.org/vocabulary/#split).*

By default, this program will give you the sequences of your splits, but will look at coverage data in terms of the parent contig. If you want to get coverage information for your splits, use `--splits-mode`. Alternatively, you can ask the program to `--report-contigs` to look at contig sequences instead. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-splits-and-coverages.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-splits-and-coverages) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
