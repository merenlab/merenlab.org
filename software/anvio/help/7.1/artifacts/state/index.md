---
layout: page
title: state [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/state
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-import-state](../../programs/anvi-import-state)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-delete-state](../../programs/anvi-delete-state)</span> <span class="artifact-r">[anvi-export-state](../../programs/anvi-export-state)</span></p>


## Description

A state describes the configuration of the anvi'o <span class="artifact-n">[interactive](/software/anvio/help/7.1/artifacts/interactive)</span> interface (i.e. the cosmetic and organizational settings that you have enabled). 

From the interface, the bottom section of the left panel enables you to save and load states. You also have the option to import states with <span class="artifact-n">[anvi-import-state](/software/anvio/help/7.1/programs/anvi-import-state)</span> or export them with <span class="artifact-n">[anvi-export-state](/software/anvio/help/7.1/programs/anvi-export-state)</span>. You can also delete states you no longer need anymore with <span class="artifact-n">[anvi-delete-state](/software/anvio/help/7.1/programs/anvi-delete-state)</span>. 

Here is the information stored in a state:
* The current item (see <span class="artifact-n">[misc-data-items](/software/anvio/help/7.1/artifacts/misc-data-items)</span>) and layers (see <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span>) displayed
    * related information, like the minimum and maximum value for the data displayed in each layer, 
* The current items order (see <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7.1/artifacts/misc-data-items-order)</span>) and layers order (see <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7.1/artifacts/misc-data-layer-orders)</span>)
* The views you have available 
* Any sample groups you have 
* Various cosmetic settings, like font size, angles, dimensions, colors, whether or not labels are displayed, etc. 
    * This includes whether your display is in circles or rectangles 
    
No more having to manually set parameters like your layer order for each bin you look at! Just save a state when the interface is adjusted to your liking, and using anvi'o will be that much easier. 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/state.md) to update this information.

