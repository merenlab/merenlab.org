---
layout: page
title: svg [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/svg
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/SVG.png" alt="SVG" style="width:100px; border:none" />

A SVG-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-display-contigs-stats](../../programs/anvi-display-contigs-stats)</span> <span class="artifact-p">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-p">[anvi-interactive](../../programs/anvi-interactive)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

SVG stands for scalable vector graphics, which is a vector-based image format. In anvi'o programs that give you pretty-looking outputs will also give you an svg, so you can look at the beauty of your data without having to open anvi'o for analysis. This also makes it easier to share with others (so you don't have to use screenshots in your poster).

As of now, <span class="artifact-n">[anvi-display-contigs-stats](/software/anvio/help/7.1/programs/anvi-display-contigs-stats)</span>, <span class="artifact-n">[anvi-display-pan](/software/anvio/help/7.1/programs/anvi-display-pan)</span>, and <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> will give you an svg output every time you click the little save button at the bottom-left corner of the settings panel in the interface. You can even customize the location of that output using the flag `--export-svg`. 

Take a look at [this blogpost](http://merenlab.org/2016/10/27/high-resolution-figures/) for an outline of how to get this svg file into a publication-quality figure. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/svg.md) to update this information.

