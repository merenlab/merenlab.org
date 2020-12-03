---
layout: page
title: interactive [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DISPLAY.png" alt="DISPLAY" style="width:100px; border:none" />

A DISPLAY-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-3dev](../../programs/anvi-3dev)</span> <span class="artifact-p">[anvi-display-contigs-stats](../../programs/anvi-display-contigs-stats)</span> <span class="artifact-p">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-p">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-p">[anvi-inspect](../../programs/anvi-inspect)</span> <span class="artifact-p">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-p">[anvi-script-snvs-to-interactive](../../programs/anvi-script-snvs-to-interactive)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This page describes the many interactive interfaces that are utilized in Anvi'o. The most well-known and sophisticated of these are the beautiful concentric circles (though they can also be displayed in other shapes) given by <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> and <span class="artifact-n">[anvi-display-pan](/software/anvio/help/main/programs/anvi-display-pan)</span>. 

If you're new to the anvi'o interactive interface, you'll probably want to check out [this tutorial for beginners](http://merenlab.org/tutorials/interactive-interface/) or the other resources on the  <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> page. 

However, there are more interfaces availible in anvi'o than just that one, so let's list them out: 

- <span class="artifact-n">[anvi-3dev](/software/anvio/help/main/programs/anvi-3dev)</span> lets you examine specific protein strcutures, along with SCV and SAAVs within it. (It even has [its own software page.](http://merenlab.org/software/anvi-3dev/). It's kind of a big deal.)

- <span class="artifact-n">[anvi-display-contigs-stats](/software/anvio/help/main/programs/anvi-display-contigs-stats)</span> shows you various stats about the contigs within a <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, such as their hmm-hits, lengths, N and L statistics, and so on.

- <span class="artifact-n">[anvi-display-metabolism](/software/anvio/help/main/programs/anvi-display-metabolism)</span> is still under development but will allow you to interactively view metabolism estimation data using <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span> under the hood. 

- <span class="artifact-n">[anvi-display-pan](/software/anvio/help/main/programs/anvi-display-pan)</span> displays information about the gene clusters that are stored in a <span class="artifact-n">[pan-db](/software/anvio/help/main/artifacts/pan-db)</span>. It lets you easily view your core and accessory genes, and can even be turned into a metapangenome through importing additional data tables. 

- <span class="artifact-n">[anvi-inspect](/software/anvio/help/main/programs/anvi-inspect)</span> lets you look at a single split across your samples, as well as the genes identified within it. This interface can also be opened from the <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> interface by asking for details about a specific split. 

- <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> displays the information in a <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span>. It lets you view the distribution of your contigs across your samples, manually bin metagenomic data into MAGSs (and refine those bins with <span class="artifact-n">[anvi-refine](/software/anvio/help/main/programs/anvi-refine)</span>), and much more. You can also use this to look at your genes instead of your contigs or [examine the genomes after a phylogenomic anlysis](http://merenlab.org/2017/06/07/phylogenomics/). Just look at that program page for a glimpse of this program's amazingness. 

- <span class="artifact-n">[anvi-script-snvs-to-interactive](/software/anvio/help/main/programs/anvi-script-snvs-to-interactive)</span> lets you view a comprehensive summary of the SNVs, SCVs, and SAAVs within your contigs. 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/interactive.md) to update this information.

