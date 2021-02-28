---
layout: page
title: "Citing anvi'o like a pro"
excerpt: "How to make sure credit goes where it is due"
categories: [anvio]
comments: true
image:
  featurefromroot: /software/anvio/help/main/images/header.png
  display: true
---

{% include _toc.html %}

Anvi'o is an evolving software ecosystem, and its components are often described in multiple studies. Thus, the best practice for your study may be to cite multiple publications if it benefits from multiple anvi'o features. 

We know that finding the best studies to cite can be a lot of work. The purpose of this page is to offer up-to-date suggestions to help you find out how to finalize your citations regarding anvi'o. But if you are unsure, please feel free to drop us a line, or find us on Slack.

{% include _join-anvio-slack.html %}

{:.notice}
Anvi'o often uses third-party software or resources (such as HMMER, Prodigal, MCL, GTDB, or NCBI) and the platform typically guides you to cite relevant work when they are used for an anvi'o analysis. Suggestions on this page are specific to anvi'o, and do not include third-party software that you should also make sure to cite properly.

We know this is difficult work and we are thankful for your attention.


## Default citation

If you have used anvi'o for anything at all please consider citing this work as it describes the software ecosystem in general which currently sits on more than 120,000 lines of code, which means any given anvi'o program benefits from the entirety of this ecosystem:

<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.1038/s41564-020-00834-3"></div>
<div class="__dimensions_badge_embed__" data-doi="10.1038/s41564-020-00834-3" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <h3><a href=" https://doi.org/10.1038/s41564-020-00834-3" target="_new">Community-led, integrated, reproducible multi-omics with anvi'o</a></h3>
    <span class="pub-authors"><span class="none">Eren AM</span>, <span class="none">Kiefl E</span>, <span class="none">Shaiber A</span>, <span class="none">Veseli I</span>, <span class="none">Miller SE</span>, <span class="none">Schechter MS</span>, <span class="none">Fink I</span>, <span class="none">Pan JN</span>, <span class="none">Yousef M</span>, <span class="none">Fogarty EC</span>, <span class="none">Trigodet F</span>, <span class="none">Watson AR</span>, <span class="none">Esen ÖC</span>, Moore RM, Clayssen Q, Lee MD, Kivenson V, Graham ED, Merrill BD, Karkman A, Blankenberg D, Eppley JM, Sjödin A, Scott JJ, Vázquez-Campos X, McKay LJ, McDaniel EA, Stevens SLR, Anderson RE, Fuessel J, Fernandez-Guerra A, Maignien L, Delmont TO, Willis AD</span>
    <span class="pub-journal"><b>Nature Microbiology</b>, 6(1):3:6.</span>
</div>

The rest of the citations on this page are specific for certain anvi'o features.


## Functional enrichment analyses

This feature was introduced by Amy Willis and Alon Shaiber, and described for the first time in this study. In addition to the default anvi'o citation, please  consider citing this work.

<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.1186/s13059-020-02195-w"></div>
<div class="__dimensions_badge_embed__" data-doi="10.1186/s13059-020-02195-w" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <h3><a href=" https://doi.org/10.1186/s13059-020-02195-w" target="_new">Functional and genetic markers of niche partitioning among enigmatic members of the human oral microbiome</a></h3>
    <span class="pub-authors"><span class="none">Shaiber A</span>, Willis AD, Delmont TO, Roux S, Chen L, <span class="none">Schmid AC</span>, <span class="none">Yousef M</span>, <span class="none">Watson AR</span>, <span class="none">Lolans K</span>, <span class="none">Esen ÖC</span>, <span class="none">Lee STM</span>, Downey N, Morrison HG, Dewhirst FE, Mark Welch JL<sup>‡</sup>, <span class="none">Eren AM<sup>‡</sup></span></span>
    <span class="pub-co-first-authors"><sup>‡</sup>Co-senior authors</span>
    <span class="pub-journal"><b>Genome Biology</b>, 21:292.</span>
</div>

In a recent study, we cited this work the following way:

> (...)
>
> **Functional enrichment analyses**. The statistical approach for enrichment analysis is defined elsewhere ([Shaiber et al. 2020](https://doi.org/10.1186/s13059-020-02195-w)), but briefly the program `anvi-compute-functional-enrichment` determined enrichment scores for functions (or metabolic modules) within groups of genomes by fitting a binomial generalized linear model (GLM) to the occurrence of each function (or complete metabolic module) in each group, and then computing a Rao test statistic, uncorrected p-values, and corrected q-values. We considered any function or metabolic module with a q-value less than 0.05 to be 'enriched' in its associated group (...)


## Snakemake workflows

There is not yet a published study that describes [our workflows](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/), although they were first introduced in the following work:

<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.1186/s13059-020-02195-w"></div>
<div class="__dimensions_badge_embed__" data-doi="10.1186/s13059-020-02195-w" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <h3><a href=" https://doi.org/10.1186/s13059-020-02195-w" target="_new">Functional and genetic markers of niche partitioning among enigmatic members of the human oral microbiome</a></h3>
    <span class="pub-authors"><span class="none">Shaiber A</span>, Willis AD, Delmont TO, Roux S, Chen L, <span class="none">Schmid AC</span>, <span class="none">Yousef M</span>, <span class="none">Watson AR</span>, <span class="none">Lolans K</span>, <span class="none">Esen ÖC</span>, <span class="none">Lee STM</span>, Downey N, Morrison HG, Dewhirst FE, Mark Welch JL<sup>‡</sup>, <span class="none">Eren AM<sup>‡</sup></span></span>
    <span class="pub-co-first-authors"><sup>‡</sup>Co-senior authors</span>
    <span class="pub-journal"><b>Genome Biology</b>, 21:292.</span>
</div>


In a recent study, we cited our workflows the following way:

> (...)
>
> **‘Omics workflows**. Whenever applicable, we automated and scaled our ‘omics analyses using the bioinformatics workflows implemented by the program `anvi-run-workflow` ([Shaiber et al. 2020](https://doi.org/10.1186/s13059-020-02195-w)) in anvi’o ([Eren et al. 2021](https://doi.org/10.1038/s41564-020-00834-3)). Anvi’o workflows implement numerous steps of bioinformatics tasks including short-read quality filtering, assembly, gene calling, functional annotation, hidden Markov model search, metagenomic read-recruitment, metagenomic binning, pangenomics, and phylogenomics. Workflows use Snakemake ([Köster and Rahmann 2012](https://doi.org/10.1093/bioinformatics/bts480)) and a tutorial is available at the URL [http://merenlab.org/anvio-workflows/](http://merenlab.org/anvio-workflows/). The following sections detail these steps.
> 
> (...)


## Metabolic reconstruction

There is not yet a published study that describes anvi'o metabolic reconstruction capabilities, although in a recent study we mentioned the program we mentioned this framework the following way:

> (...)
>
> **Analysis of metabolic modules and enrichment**. We calculated the level of completeness for a given KEGG module ([Kanehisa et al. 2014](https://doi.org/10.1093/nar/gkt1076); [Kanehisa et al. 2017](https://doi.org/10.1093/nar/gkw1092)) in our genomes using the program `anvi-estimate-metabolism`, which leveraged previous annotation of genes with KEGG orthologs (KOs) (see the section ‘Processing of contigs’). Then, the program `anvi-compute-functional-enrichment` determined whether a given metabolic module was enriched in based on the output from `anvi-estimate-metabolism`.  The URL [https://merenlab.org/m/anvi-estimate-metabolism](https://merenlab.org/m/anvi-estimate-metabolism) serves a tutorial for this program which details the modes of usage and output file formats (...)


## Single-amino acid variants

If you are using anvi'o to study microbial population genetics through single-codon variants or single-amino acid variants, please consider also citing this work:

<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.7554/eLife.46497"></div>
<div class="__dimensions_badge_embed__" data-doi="10.7554/eLife.46497" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <h3><a href=" https://doi.org/10.7554/eLife.46497" target="_new">Single-amino acid variants reveal evolutionary processes that shape the biogeography of a global SAR11 subclade</a></h3>
    <span class="pub-authors"><span class="none">Delmont TO<sup>☯</sup></span>, <span class="none">Kiefl E<sup>☯</sup></span>, Kilinc O, <span class="none">Esen ÖC</span>, Uysal I, Rappé MS, Giovannoni S, <span class="none">Eren AM</span></span>
    <span class="pub-co-first-authors"><sup>☯</sup>Co-first authors</span>
    <span class="pub-journal"><b>eLife</b>, 8:e46497.</span>
</div>


## Metapangenomics

The metapangenomics was first introduced in this study. If you are using anvi'o to investigate how to bring together pangenomes and metagenomes, please consider citing this work as well.

<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.7717/peerj.4320"></div>
<div class="__dimensions_badge_embed__" data-doi="10.7717/peerj.4320" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <h3><a href=" https://doi.org/10.7717/peerj.4320" target="_new">Linking pangenomes and metagenomes: the Prochlorococcus metapangenome</a></h3>
    <span class="pub-authors"><span class="none">Delmont TO</span>, <span class="none">Eren AM</span></span>
<span class="pub-journal"><b>PeerJ</b>, 6:e4320.</span>
</div>


## Metagenomic binning / refinement only

If you used anvi'o only for metagenomic binning or for the refinement of genomes, please consider citing this study, too.

<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.7717/peerj.1319"></div>
<div class="__dimensions_badge_embed__" data-doi="10.7717/peerj.1319" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <h3><a href=" https://doi.org/10.7717/peerj.1319" target="_new">Anvi’o: an advanced analysis and visualization platform for ‘omics data</a></h3>
    <span class="pub-authors"><span class="none">Eren AM</span>, <span class="none">Esen ÖC</span>, Quince C, Vineis JH, Morrison HG, Sogin ML, <span class="none">Delmont TO</span></span>
    <span class="pub-journal"><b>PeerJ</b>, 3:e1319.</span>
</div>