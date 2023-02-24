---
layout: project-anvio-structure-main
title: "Integrating metagenomic sequence variants and predicted protein structures"
redirect_from:
  - /anvio-structure
  - /projects/anvio-structure/
excerpt: "Anvi'o Structure Project Page"
image:
  feature: eel-pond.jpg
  display: true
---

{% include _toc.html %}

# Anvi'o structure: Integrating metagenomic sequence variants and predicted protein structures

Anvi'o structure is an integrated component of anvi'o focused on the analysis and visualization of
single-amino acid and single-codon variants derived from metagenomic/metatranscriptomic data in the
context of predicted protein structures and binding sites. Its purpose is to enable microbial
ecologists with the tools necessary to investigate microbial diversity within the context of the
structural properties of translated gene products, thereby bridging a gap between microbial ecology
and structural biology.

# Install

Anvi'o structure is an integrated part of the anvi'o framework. By **[installing anvi'o](/install-anvio/)**,
you are installing everything that makes up anvi'o structure, as well.

# How to cite

If you use anvi'o structure, please cite this paper:


<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.1126/sciadv.abq4632"></div>
<div class="__dimensions_badge_embed__" data-doi="10.1126/sciadv.abq4632" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <span class="pub-title"><a href="https://doi.org/10.1126/sciadv.abq4632" target="_new">Structure-informed microbial population genetics elucidate selective pressures that shape protein evolution</a></span>
    <span class="pub-authors"><span class="pub-member-author">Kiefl E</span>, Esen ÖC, <span class="pub-member-author">Miller SE</span>, Kroll KL, Willis AD, Rappé MS, Pan T, <span class="pub-member-author">Eren AM</span></span>
    <div class="pub-info">
    <div class="pub-featured-image">
    <a href="/images/pubs/kiefl_et_al_metagenomics_plus_protein_structures.png"><img src="/images/pubs/kiefl_et_al_metagenomics_plus_protein_structures.png" style="max-width: 100px; max-height: 80px; width: auto; border: none; height: auto; margin: 0 auto; display: block; transform: translateY(15%);" /></a>
    </div>
    <div class="pub-highlights">
    <span style="display: inline-block; padding-bottom: 5px;">- A study that describes an approach to integrate <b>environmental microbiology</b> with recent advances in <b>protein structure prediction</b>, and illustrates the tight association between <b>intra-population genetic variants</b>, <b>environmental selective pressures</b>, and <b>structural properties of proteins</b>.</span><br /><span style="display: inline-block; padding-bottom: 5px;">- Demonstrates a quantifiable link between (1) the magnitude of selective pressures over key metabolic <b>genes</b> (e.g., glutamine synthase of the central nitrogen metabolism), (2) the availability of key <b>nutrients</b> in the environment (e.g., nitrate), and (3) the maintenance of nonsynonymous <b>variants</b> near protein active sites.</span><br /><span style="display: inline-block; padding-bottom: 5px;">- Shows that the interplay between selective pressures and protein structures also maintains <b>synonymous variants</b> -- revealing a quantifiable link between <b>translational accuracy</b> and fluctuating <b>selective pressures</b>.</span><br /><span style="display: inline-block; padding-bottom: 5px;">- Comes with a <a href="https://merenlab.org/data/anvio-structure/chapter-I/"><b>reproducible bioinformatics workflow</b></a> that offers detailed access to computational steps used in the study that spans from metagenomic read recruitment and profiling to the integration of environmental variants and predicted protein structures.</span>
    </div>
    </div>
    <span class="pub-journal"> 📚 <b>Science Advances</b>, 9(8):eabq4632 | 🔍 <a href="http://scholar.google.com/scholar?hl=en&amp;q=Structure-informed+microbial+population+genetics+elucidate+selective+pressures+that+shape+protein+evolution" target="_blank">Google Scholar</a> | 🔗 <a href="https://doi.org/10.1126/sciadv.abq4632" target="_blank">doi:10.1126/sciadv.abq4632</a></span>
</div>

# Made possible by

Like all software, anvi'o structure is built upon other software, all of which is open source. It is
especially indebted to the following list of projects, without which anvi'o structure would not have been
developed:

### MODELLER

{:.notice}
**Citation**: [doi:10.1002/0471250953.bi0506s15](https://doi.org/10.1002/0471250953.bi0506s15)

{:.notice}
**Citation**: [doi:10.1146/annurev.biophys.29.1.291](https://doi.org/10.1146/annurev.biophys.29.1.291)

{:.notice}
**Citation**: [doi:10.1006/jmbi.1993.1626](https://doi.org/10.1006/jmbi.1993.1626)

MODELLER is the program anvi'o uses to predict protein structure based on experimentally solved structures in the Protein Data Bank. We'll talk more specifically about how it accomplishes this in the following section, but for now you need to make sure it's installed on your computer. For that, check out these instructions to see if you have it installed ([click me](http://merenlab.org/2016/06/18/installing-third-party-software/#modeller)), and how to install it if you don't. We've tried to make it as simple for you as possible.

### NGL

{:.notice}
**Citation**: [doi:10.1093/nar/gkv402](https://doi.org/10.1093/nar/gkv402)

{:.notice}
**Citation**: [doi:10.1093/bioinformatics/bty419](https://doi.org/10.1093/bioinformatics/bty419)

NGL (NGL) is an open-source project for visualizing biomolecules. This browser-based solution to visualization means you don't have to install anything, and you can thank them for that. Özcan has been the mastermind behind incorporating NGL's visualizations into anvi'o, and says he continues to be impressed with their excellent code, documentation, and open-source approach to science.


# Posts and tutorials about anvi'o structure
