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

If you use anvi'o structure, please cite this paper.

Kiefl E, Esen ÖC, Eren AM (2021) FIXME. (*in preparation*)

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
