---
layout: page
title: interactive [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/interactive
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DISPLAY.png" alt="DISPLAY" style="width:100px; border:none" />

A DISPLAY-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-display-contigs-stats](../../programs/anvi-display-contigs-stats)</span> <span class="artifact-p">[anvi-display-functions](../../programs/anvi-display-functions)</span> <span class="artifact-p">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-p">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-p">[anvi-display-structure](../../programs/anvi-display-structure)</span> <span class="artifact-p">[anvi-inspect](../../programs/anvi-inspect)</span> <span class="artifact-p">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-p">[anvi-script-checkm-tree-to-interactive](../../programs/anvi-script-checkm-tree-to-interactive)</span> <span class="artifact-p">[anvi-script-gen-functions-per-group-stats-output](../../programs/anvi-script-gen-functions-per-group-stats-output)</span> <span class="artifact-p">[anvi-script-snvs-to-interactive](../../programs/anvi-script-snvs-to-interactive)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This page describes general properties of anvi'o interactive displays and programs that offer anvi'o interactive artifacts.

## Terminology

Anvi'o uses a simple terminology to address various aspects of interactive displays it produces, such as items, layers, views, orders, and so on. The purpose of this section is to provide some insights into these terminology using the figure below:

![an anvi'o display](../../images/interactive_interface/anvio_display_template.png){:.center-img}

Even though the figure is a product of <span class="artifact-n">[anvi-display-pan](/software/anvio/help/7.1/programs/anvi-display-pan)</span>, the general terminology does not change across different interfaces, including the default visualizations of <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span>. Here are the descriptions of numbered areas in the figure:

* The tree denoted by **(1)** shows the organization of each `item`. Items could be contigs, gene clusters, bins, genes, or anything else depending on which mode the anvi'o interactive interface was initiated. The structure that orders items and denoted by **(1)** in the figure can be a phylogenetic or phylogenomic tree, or a dendrogram produced by a hierarchical clustering algorithm. In addition, there may be nothing there, if the user has requested or set a linear items order through <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7.1/artifacts/misc-data-items-order)</span>.
* Each concentric circle underneath the number **(2)** is called a `layer` and the data shown for items and layers as a whole is called a `view`. A **layer** can be a genome, a metagenome, or anything else depending on which mode the anvi'o interactive was initiated. The **view** is like a data table where a datum is set for each **item** in each **layer**. The view data is typically computed by anvi’o and stored in pan databases by <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/7.1/programs/anvi-pan-genome)</span> or profile databases by <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span>. The user add another view to the relevant combo box in the interface by providing a TAB-delimited file to <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> through the command line argument `--additional-view`, or add new layers to extend these vies with additional data through <span class="artifact-n">[misc-data-items](/software/anvio/help/7.1/artifacts/misc-data-items)</span>.
* The tree denoted by **(3)** shows a specific ordering of layers. Anvi'o will compute various layer orders automatically based on available **view** depending on the analysis or visualization mode, and users can extend available **layer orders** through <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/7.1/artifacts/misc-data-layer-orders)</span>.
* What is shown by **(4)** is the additional data for layers. the user can extend this section with additional information on layers using the <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span>.

The orchestrated use of <span class="artifact-n">[anvi-import-misc-data](/software/anvio/help/7.1/programs/anvi-import-misc-data)</span>, <span class="artifact-n">[anvi-export-misc-data](/software/anvio/help/7.1/programs/anvi-export-misc-data)</span>, and <span class="artifact-n">[anvi-delete-misc-data](/software/anvio/help/7.1/programs/anvi-delete-misc-data)</span> provides a powerful framework to decorate items or layers in a display and enhance visualization of complex data. Please take a look at the following article on how to extend anvi'o displays:

* [https://merenlab.org/2017/12/11/additional-data-tables/](https://merenlab.org/2017/12/11/additional-data-tables/)

## Programs that give interactive access

If you're new to the anvi'o interactive interface, you'll probably want to check out [this tutorial for beginners](http://merenlab.org/tutorials/interactive-interface/) or the other resources on the  <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> page.

However, there are more interfaces available in anvi'o than just that one, so let's list them out:

- <span class="artifact-n">[anvi-display-structure](/software/anvio/help/7.1/programs/anvi-display-structure)</span> lets you examine specific protein structures, along with SCV and SAAVs within it. (It even has [its own software page.](http://merenlab.org/software/anvio-structure/). It's kind of a big deal.)

- <span class="artifact-n">[anvi-display-contigs-stats](/software/anvio/help/7.1/programs/anvi-display-contigs-stats)</span> shows you various stats about the contigs within a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, such as their hmm-hits, lengths, N and L statistics, and so on.

- <span class="artifact-n">[anvi-display-functions](/software/anvio/help/7.1/programs/anvi-display-functions)</span> lets you quickly browse the functional pool for a given set of genomes or metagenomes.

- <span class="artifact-n">[anvi-display-metabolism](/software/anvio/help/7.1/programs/anvi-display-metabolism)</span> is still under development but will allow you to interactively view metabolism estimation data using <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7.1/programs/anvi-estimate-metabolism)</span> under the hood.

- <span class="artifact-n">[anvi-display-pan](/software/anvio/help/7.1/programs/anvi-display-pan)</span> displays information about the gene clusters that are stored in a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>. It lets you easily view your core and accessory genes, and can even be turned into a metapangenome through importing additional data tables.

- <span class="artifact-n">[anvi-inspect](/software/anvio/help/7.1/programs/anvi-inspect)</span> lets you look at a single split across your samples, as well as the genes identified within it. This interface can also be opened from the <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> interface by asking for details about a specific split.

- <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> displays the information in a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>. It lets you view the distribution of your contigs across your samples, manually bin metagenomic data into MAGSs (and refine those bins with <span class="artifact-n">[anvi-refine](/software/anvio/help/7.1/programs/anvi-refine)</span>), and much more. You can also use this to look at your genes instead of your contigs or [examine the genomes after a phylogenomic analysis](http://merenlab.org/2017/06/07/phylogenomics/). Just look at that program page for a glimpse of this program's amazingness.

- <span class="artifact-n">[anvi-script-snvs-to-interactive](/software/anvio/help/7.1/programs/anvi-script-snvs-to-interactive)</span> lets you view a comprehensive summary of the SNVs, SCVs, and SAAVs within your contigs.

## Artifacts that give interactive access

- <span class="artifact-n">[gene-cluster-inspection](/software/anvio/help/7.1/artifacts/gene-cluster-inspection)</span> lets you examine specific gene clusters.

- <span class="artifact-n">[contig-inspection](/software/anvio/help/7.1/artifacts/contig-inspection)</span> shows you detailed contig information.

## An overview of the display

The interactive interface has three major areas of interaction:

* The space for visualization in the middle area,
* The Settings panel on the left of the screen,
* And three additional panels for 'News', 'Description', and 'Mouse' on the right.

Each panel is important, but the most important and functionally rich one is the 'Settings' panel.

### Settings panel

If closed, the settings panel can be opened by clicking on the little button on the left-middle part of your browser. When opened, you will see multiple tabs:

![an anvi'o settings panel](../../images/interactive_interface/interactive-settings-panel-tabs.png){:.center-img}

But before we start talking about these tabs, it is worthwhile to mention that at the bottom of the settings panel you will find a section with tiny controls that are available in all tabs:

![settings panel bottom controls](../../images/interactive_interface/interactive-settings-bottom.png){:.center-img}

Through these controls you can,

* __Create or refresh__ the display when necessary using the draw button (some changes require you to do that),

* __Zoom in, zoom out, and center__ the display.

* __Download your display as an SVG file.__

Finally, at the top-right of the Settings panel header you will find a dropdown menu (hamburger menu) which provides links to external information, resources, and issue-reporting related to anvi'o.

OK. Let's talk about each tab you will find in the settings panel.


### Main Tab

This is one of the most frequently used tabs in the interface, and there are multiple sections in it (keeps growing over time, so things may be missing here).

![an anvi'o main tab](../../images/interactive_interface/interactive-settings-display-additional-settings.png){:.center-img}

* **Display subsection**. Provides high level options for adjusting _items order_, _view_, and _drawing type_.

Clicking the _Show Additional Settings_ button provides access to myriad additional, more-granular adjustments, including,

* **Dendrogram subsection**. _Radius_ and _Angle_ , and _Edge length normalization_ adjustments for the dendrogram.
* **Branch support subsection**. Settings for displaying _bootstrap values_ on the dendrogram.
* **Selections subsection**. To adjust _height_, _grid_ and/or _shade_ display, as well as selection _name_ settings.
* **Layers subsection**. Display and label settings.
* **Performance subsection**. Whether the SVG output is optimized for performance or granularity (very advanced stuff).
* **Layers subsection**. This is arguably the most important subsection in the Main tab that enables you to make very precise adjustments to how things should look like on your screen. You can adjust individual layer attributes like _color_, display _type_, _height_  and _min/max_ values. Click + drag each layer to rearrange how layers are ordered. Or _edit attributes for multiple layers_ as well.

![an anvi'o settings layers](../../images/interactive_interface/interactive-settings-layers.png){:.center-img}

Mastering these in the Main Tab will minimize the post-processing of your anvi'o figures for high-quality and good-looking publication ready images.

### Layers tab

Through the layers tab you can,

- __Change general settings for the tree__ (i.e., switching between circle or rectengular displays, changing tree radius or width), __and layers__ (i.e., editing layer margins, or activating custom layer margins).

- __Load or save states__ to store all visual settings, or load a previously saved state.

- __Customize individual__ layers by switching between different __display modes__ depending on the layer type (i.e., ‘text’ or ‘color’ mode for categorical layers, or ‘bar’ or ‘intensity’ mode for numerical layers), __set normalization__ (i.e., ‘square-root’, or ‘log’ normalization), __minimum, and maximum cutoff__ values for numerical layers, or set __layer height__, and __layer margin__ (i.e., its distance from the previous layer).

- Use the __multi-selector__ at the bottom to change settings for multiple layers at once.

![an anvi'o layers tab](../../images/interactive_interface/interactive-settings-layers-tab.png){:.center-img}


### Samples tab

Samples tab is for the additional data you provide the interface through a samples database (see samples order and samples information sections above). Through this layer you can,

- __Change the order__ of layers using automatically-generated or user-provided orders of layers using the Sample order combo box,

- __Customize individual samples information entries.__ Changes in this tab can be reflected to the current display without re-drawing the entire tree unless the sample order is changed.

### Bins tab

Anvi’o allows you to create selections of items shown in the display (whether they are contigs, gene clusters, or any other type of data shown in the display). Bins tab allow you to maintain these selections. Any selection on the tree will be added to active bin in this tab (the state radio button next to a bin defines its activity). Through this tab you can,

- __Create or delete bins, set bin names, change the color of a given bin__, or sort bins based on their name, the number of units they carry, or completion and contamination estimates (completion / contamination estimates are only computed for genomic or metagenomic analyses).

- View the __number of selected units__ in a given bin, and see the __list of names in the selection__ by clicking the button that shows the number of units described in the bin.

- __Store a collection of bins__, or __load a previously stored collection.__

![an anvi'o bins tab](../../images/interactive_interface/interactive-settings-bins-tab.png){:.center-img}

### Legends tab

The legends tab enables users to easily change individual or batch legend colors for any of their additional data items
<!-- grab legend example from infant gut w/ additional data layer -->

![an anvi'o legends tab](../../images/interactive_interface/interactive-settings-legends-tab.png){:.center-img}

### Search tab

It does what the name suggests. Using this tab you can,

- __Build expressions to search items__ visualized in the main display.

- __Highlight matches__, and __append__ them to, or __remove__ them from the __selected bin__ in the Bins tab.



### Mouse panel

The mouse panel displays the value of items underneath the mouse pointer while the user browse the tree.

Displaying the numerical or categorical value of an item shown on the tree is not an easy task. We originally thought that displaying pop-up windows would solve it, but besides the great overhead, it often became a nuisance while browsing parts of the tree. We could show those pop-up displays only when use clicks on the tree, however click-behavior is much more appropriate to add or remove individual items from a bin, hence, it wasn’t the best solution either. So we came up with the ‘mouse panel’. You have a better idea? I am not surprised! We would love to try improve your experience: please enter an issue, and let’s discuss.

### News panel
The news panel provides information and external links tracking major Anvi'o releases and development updates.

### Description panel

- The description panel is a flexible, multipurpose space where users can,
- Store notes, comments, and any other stray items related to their project, in a feature-rich markdown environment.
- Display context, references, reproducibility instructions, and any other salient details for published figures.
![The Description panel in action](../../images/interactive_interface/interactive-settings-description-panel.png){:.center-img}

## Interactive interface tips + tricks

Here are some small conveniences that may help the interface serve you better (we are happy to expand these little tricks with your suggestions).

* You can zoom to a section of the display by making a rectangular selection of the area __while the pressing the shift button.__

* You can click an entire branch to add items into the selected bin, and remove them by __right-clicking__ a branch.

* If you click a branch __while pressing the `Command` or `CTRL` button__, it will create a new bin, and add the content of the selection into that bin.

* Tired of selecting items for binning one by one? __right-click__ on an item and select __Mark item as 'range start'__ to set an 'in point', then __right-click__ on another item and select __Add items in range to active bin__ or __Remove items in range from any bin__ to manipulate many items with few clicks. Nice!

* By pressing `1`,`2`,`3`,`4`, and`5`, you can go between Layers, Bins, Samples, Mouse, and Search tabs!

## Keyboard shortcuts

The interactive interface recognizes a handful of keyboard shortcuts to help speed up your workflow

- The `S` key toggles the Settings panel
- The `M` key toggles the Mouse panel
- The `N` key toggles the Description panel
- The `W` key toggles the News panel
- The `D` key triggers a redraw of your visualization
- The `T` key toggles showing the Title panel
- Keys `1` through `5` will toggle between tabs within the Settings panel, granted the Settings panel is currently shown.
- `CTRL`+`Z` and `CTRL`+`SHIFT`+`Z` will undo or redo bin actions, respectively.







{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/interactive.md) to update this information.

