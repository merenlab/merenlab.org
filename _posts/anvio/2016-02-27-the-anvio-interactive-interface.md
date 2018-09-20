---
layout: post
title: "The anvi'o interactive interface"
excerpt: "Data types, usage tips, and other stuff about the interface"
modified: 2016-02-27
tags: [interactive]
categories: [anvio]
comments: true
authors: [meren, ozcan]
---

{% include _toc.html %}

<p style="text-align: right; font-style: italic; color: #AAAAAA;">This document last updated on April 1<sup>st</sup>, 2016.</p>

Anvi'o is a comprehensive 'omics platform with a large codebase to perform [wide range of tasks](https://github.com/meren/anvio/tree/master/bin){:target="_blank"}, and [in-depth analyses]({% post_url anvio/2015-12-09-musings-over-commamox %}){:target="_blank"} of large datasets.

With anvi'o you can do [metagenomic binning](https://peerj.com/articles/1319/){:target="_blank"}, characterize [single-nucleotide variation]({% post_url anvio/2015-07-20-analyzing-variability %}){:target="_blank"}, study [bacterial pangenomes]({% post_url anvio/2016-11-08-pangenomics-v2 %}){:target="_blank"}, [benchmark]({% post_url anvio/2015-06-23-comparing-different-mapping-software %}){:target="_blank"} software tools, predict [number of bacterial genomes in a metagenomic assembly](http://merenlab.org/2015/12/07/predicting-number-of-genomes/){:target="_blank"}, or even [remove contamination from eukaryotic assembly projects](https://peerj.com/preprints/1695/){:target="_blank"}. Anvi'o's 'versatility' partly comes from its integrated visualization framework that allows the user to *see* all these different types of data, and *interact* with them.

<div class="centerimg">
<img src="http://i.imgur.com/d1c7bUY.png?1" style="margin: 3px;" />
<img src="http://i.imgur.com/tIG2ZMJ.png?1" style="margin: 3px;" />
<img src="http://i.imgur.com/ILhiAbP.png?1" style="margin: 3px;" />
<img src="http://i.imgur.com/T36nS6D.png?1" style="margin: 3px;" />
<img src="http://i.imgur.com/iGuCRnu.jpg?1" style="margin: 3px;" />
</div>

The anvi'o interactive interface is a fully customizable visualization environment that is accessible through an intuitive interface to efficiently visualize complex data. It can handle large datasets, and it's [source code](http://github.com/meren/anvio){:target="_blank"} is freely available within the anvi'o platform.

Although it is fully integrated with core anvi'o operations detailed in [the metagenomic workflow tutorial]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}){:target="_blank"}, the visualization environment can be initiated in an *ad hoc* manner by using the `anvi-interactive` program with `--manual-mode` flag, or through anvi'server, without an anvi'o installation. In summary, if you have a matrix file, anvi'o may be useful to you to generate high-quality, publication-ready figures with mouse clicks.

**The purpose** of this article is to provide a more detailed description of the interface by demonstrating the data types the interface can work with, and later details of the user interface. 

---

<h3>A little note on our ongoing project, anvi'server</h3>

To make the anvi'o interactive interface more accessible, we teamed up with [Tobias Paczian](https://github.com/paczian){:target="_blank"}, and with his remarkable efforts created a web service. This new service, which we call **anvi'server**, is now running at [http://anvi-server.org](http://anvi-server.org){:target="_blank"}. Through anvi'server, you can perform anvi'o visualizations by uploading your data through a simple interface:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/upload-project.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/upload-project.png" width="30%" /></a>
</div>

Topics covered for the remainder of this article are directly applicable to the interactive interface whether it is accessed through local anvi'o installations, or through the anvi'server.

{:.notice}
Please note that the anvi'server is under active development, and your testing efforts will be greatly appreciated. Please don't hesitate to get in touch with us if you have any questions.

---

# 1. Data types and input files

 
The purpose of this section is to provide examples for each of the data type (and input files) the interactive interface can work with. For this, I will start with a simple tree, and add layers step by step to describe different data types.

You can follow these examples in two ways:

* **Using anvi-interactive in your terminal**: For each data type I will either provide a link to the files used in the example command line, or give an example file structure so you can try them on your own files.
* **Using [http://anvi-server.org](http://anvi-server.org){:target="_blank"}**: The other option is to use our new anvi'server without installing anvi'o. If your only purpose with the interactive interface is to do an *ad hoc* visualization, I think this would be the best way to go. Otherwise you can read about the [ways to install the platform]({% post_url anvio/2016-06-26-installation-v2 %}) on your own server or laptop.

{:.notice}
Command lines mentioned in this article are run on anvi'o version 2 or later. You can check your verison using `anvi-profile -v`.

OK. Let's start.

## 1.1 Newick tree

The least you can do with the anvi'o interactive interface would be to visualize a newick-formatted tree. The [tree file](https://github.com/meren/anvio/blob/master/tests/sandbox/files_for_manual_interactive/tree.txt){:target="_blank"} I use for this example contains 300 leaves. This is how I run the interactive interface from the command line:

{% highlight bash %}
anvi-interactive -t tree.txt -p profile.db --title 'Interface Demo I' --manual
{% endhighlight %}

Which gives me this:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/01.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/01.png" width="50%" /></a>
</div>

{:.notice}
Here is the anvi'server link for this visualization: [http://anvi-server.org/public/meren/interface_demo_I](http://anvi-server.org/public/meren/interface_demo_I){:target="_blank"}


## 1.2 Numerical

Let's assume you for each item you have in the previous tree, you have multiple numerical values you want to overlay on the tree in a TAB-delimited matrix file that looks [like this](https://github.com/meren/anvio/blob/master/tests/sandbox/files_for_manual_interactive/view_data.txt){:target="_blank"}:

|contig|c1|c2|c3|
|:--|:--:|:--:|:--:|
|cathetus|13.66596066|9.590942918|46.01477372|
|centrist|11.32571669|10.08709331|36.11828559|
|cascaded|10.82858312|6.312884813|34.50972839|
|crocking|12.46382532|6.595503878|42.43493547|
|cinchona|12.78031823|8.77463282|30.56242617|
|couchant|11.90769144|4.366432391|46.86065633|
|creviced|12.88691388|9.368377696|42.06930372|
|clovered|10.51030592|8.200407215|27.51057477|
|contempt|11.25596871|7.310381191|28.93509622|
|(...)|(...)|(...)|(...)|

This data can be visualized along with the tree this way:

{% highlight bash %}
anvi-interactive -t tree.txt -d view_data.txt -p profile.db --title 'Interface Demo II' --manual
{% endhighlight %}

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/02.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/02.png" width="50%" /></a>
</div>

{:.notice}
Here is the anvi'server link for this visualization: [http://anvi-server.org/public/meren/interface_demo_II](http://anvi-server.org/public/meren/interface_demo_II){:target="_blank"}

*If you only have this matrix file but not the tree file, you can run this command to get the tree file:*

{% highlight bash %}
anvi-matrix-to-newick view_data.txt -o tree.txt
{% endhighlight %}

## 1.3 Text

Adding a text layer is as simple as adding a column of text values in your matrix file:

|contig|c1|c2|c3|text|
|:--|:--:|:--:|:--:|:--|
|backrest|11.68762604|31.81211217|16.14890468|nmwje|
|backward|8.383113248|36.27705135|12.26495265|bqmyujrpsrddoefhi|
|backwind|14.30588649|33.19818058|13.90379515|hkferlchpmzix|
|backyard|12.98490431|35.1528258|14.2861336|advoebfkyhmg|
|bacteria|6.655636411|34.45757002|13.67026608|lqmcwnhywco|
|bacterin|7.664397508|35.28860294|16.85214201|vxqdmn|
|baetylus|13.72179583|32.81162715|14.83756192|fkgpydiowgyhfxxwlpj|
|bagpiped|9.095058043|33.0406183|13.48649074|ijmnur|
|balconet|10.94823472|33.03846226|16.87202682|ecizgs|
|(...)|(...)|(...)|(...)|(...)|

Now I can run the same command on this updated matrix file:

{% highlight bash %}
anvi-interactive -t tree.txt -d view_data.txt -p profile.db --title 'Interface Demo III' --manual
{% endhighlight %}

And here is the result:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/03.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/03.png" width="50%" /></a>
</div>

{:.notice}
Anvi'server link: [http://anvi-server.org/public/meren/interface_demo_III](http://anvi-server.org/public/meren/interface_demo_III){:target="_blank"}

Please note that when I first visualized only the tree file, there were labels for each leaf. However, when I visualized the tree file along with the data, labels disappeared. Anvi'o is designed to visualize trees with thousands of items, where having labels rarely is useful. Therefore the default behavior of the visualization interface is to omit labels when there is data, and labels are shown only when the user wants to visualize a single tree. If you would like to show labels, you can simply use the text data type to duplicate the first column in your matrix file:

|contig|labels|c1|c2|c3|text|
|:--|:--|:--:|:--:|:--:|:--|
|backrest|backrest|11.68762604|31.81211217|16.14890468|nmwje|
|backward|backward|8.383113248|36.27705135|12.26495265|bqmyujrpsrddoefhi|
|backwind|backwind|14.30588649|33.19818058|13.90379515|hkferlchpmzix|
|backyard|backyard|12.98490431|35.1528258|14.2861336|advoebfkyhmg|
|bacteria|bacteria|6.655636411|34.45757002|13.67026608|lqmcwnhywco|
|bacterin|bacterin|7.664397508|35.28860294|16.85214201|vxqdmn|
|baetylus|baetylus|13.72179583|32.81162715|14.83756192|fkgpydiowgyhfxxwlpj|
|bagpiped|bagpiped|9.095058043|33.0406183|13.48649074|ijmnur|
|balconet|balconet|10.94823472|33.03846226|16.87202682|ecizgs|
|(...)|(...)|(...)|(...)|(...)|(...)|


## 1.4 Categorical

Here I am adding a column of categorical data in the same file:

|contig|c1|c2|c3|categorical|text|
|:--|:--:|:--:|:--:|:--:|:--|
|backrest|11.68762604|31.81211217|16.14890468|a|nmwje|
|backward|8.383113248|36.27705135|12.26495265|a|bqmyujrpsrddoefhi|
|backwind|14.30588649|33.19818058|13.90379515|b|hkferlchpmzix|
|backyard|12.98490431|35.1528258|14.2861336|b|advoebfkyhmg|
|bacteria|6.655636411|34.45757002|13.67026608|b|lqmcwnhywco|
|bacterin|7.664397508|35.28860294|16.85214201|c|vxqdmn|
|baetylus|13.72179583|32.81162715|14.83756192|c|fkgpydiowgyhfxxwlpj|
|bagpiped|9.095058043|33.0406183|13.48649074|c|ijmnur|
|balconet|10.94823472|33.03846226|16.87202682|c|ecizgs|
|(...)|(...)|(...)|(...)|(...)|(...)|

The same command line for the matrix file above:

{% highlight bash %}
anvi-interactive -t tree.txt -d view_data.txt -p profile.db --title 'Interface Demo IV' --manual
{% endhighlight %}

Produces this:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/04.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/04.png" width="50%" /></a>
</div>

{:.notice}
Anvi'server link: [http://anvi-server.org/public/meren/interface_demo_IV](http://anvi-server.org/public/meren/interface_demo_IV){:target="_blank"}

Some of you are probably asking themselves '*what is the difference between text and categorical data?*'. Very good question! Essentially they are both the same. If the unique number of items in a given dataset of non-integer values more than 12, the interactive interface assumes that this is a text layer, and instead of assigning random colors to each item, it shows it as such. In contrast, if there are 12 or less unique values, the interface by default treats it as a categorical data layer, and assigns random colors. Of course these colors can be changed quite easily through the interface, or programmatically by processing the 'state' file (more on this later).

## 1.5 Stacked bars

Stacked bars are a bit tricky compared to the other data types, but nothing too complicated. Here is an example addition to our data file:

|contig|c1|c2|c3|categorical|bars!A|bars!B|bars!C|text|
|:--|:--:|:--:|:--:|:--:|:--:|:--|
|backrest|11.68762604|31.81211217|16.14890468|b|278|23|1|nmwje|
|backward|8.383113248|36.27705135|12.26495265|b|249|52|2|bqmyujrpsrddoefhi|
|backwind|14.30588649|33.19818058|13.90379515|b|269|32|3|hkferlchpmzix|
|backyard|12.98490431|35.1528258|14.2861336|b|205|96|4|advoebfkyhmg|
|bacteria|6.655636411|34.45757002|13.67026608|b|263|38|5|lqmcwnhywco|
|bacterin|7.664397508|35.28860294|16.85214201|b|298|3|6|vxqdmn|
|baetylus|13.72179583|32.81162715|14.83756192|b|219|82|7|fkgpydiowgyhfxxwlpj|
|bagpiped|9.095058043|33.0406183|13.48649074|b|212|89|8|ijmnur|
|balconet|10.94823472|33.03846226|16.87202682|b|289|12|9|ecizgs|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

As you can see, the header field contains multipart information. Columns that sharing same label before exclamation mark will be grouped and will appear as a single stacked bar.

The command line for this matrix:

{% highlight bash %}
anvi-interactive -t tree.txt -d view_data.txt -p profile.db --title 'Interface Demo V' --manual
{% endhighlight %}

Produces this:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/05.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/05.png" width="50%" /></a>
</div>

And you will be able to view the values of individual bars in mouse panel with their assigned color when you hover an item.

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/08.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/08.png" /></a>
</div>


{:.notice}
Anvi'server link: [http://anvi-server.org/public/meren/interface_demo_V](http://anvi-server.org/public/meren/interface_demo_V){:target="_blank"}

## 1.6 Samples information

The *samples information* file contains one or more data types (numerical, categorical, or stacked-bar data) to provide **contextual data for layers of interest**. *Layers of interest* are, in most cases, also correspond to samples in the study. In our test case, these layers correspond to `c1`, `c2` and `c3` (see the previous matrix file).

Here is an example *samples information* file for the view data matrix file we have been using in previous steps for visualization:

|samples|numerical_01|numerical_02|categorical|stacked_bar!X|stacked_bar!Y|stacked_bar!Z|
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
|c1|100|5|A|1|2|3|
|c2|200|4|B|2|3|1|
|c3|300|3|B|3|1|2|

Please note that this time the first column is composed of layer names that appeared as rows in the data matrix file. If you are using the interactive interface via the command line, you first need to crate a samples database using this file, if you are using anvi'server, you can simply provide the TAB-delimited file via the new project window. Here is the command line (assuming the view data is identical to the one used in the previous example):

{% highlight bash %}
anvi-gen-samples-info-database --samples-information samples-information.txt -o samples.db
anvi-interactive -t tree.txt -d view_data.txt -p new_profile.db --title 'Interface Demo VI' --manual -s samples.db
{% endhighlight %}

Which produces this:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/06.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/06.png" width="50%" /></a>
</div>

{:.notice}
Anvi'server link: [http://anvi-server.org/public/meren/interface_demo_VI](http://anvi-server.org/public/meren/interface_demo_VI){:target="_blank"}

## 1.7 Samples order

The *samples order* file contains information about different *orderings* of samples, so the user can access to this information from the interface to arrange **layers of interest**. An order can be a comma separated list of sample names, or it can be a newick tree for the organization of samples.

A *samples order* file can be used with or without a *samples information* file (see `anvi-gen-samples-db -h` for help if you are following these examles from your terminal). Here is an example *samples order* file for the data matrix file we have been using:

|attributes|basic|newick|
|:--|:--:|:--:|
|test_tree||(c2:0.0370199,(c1:0.0227268,c3:0.0227268)Int3:0.0370199);|
|test_list|c3,c1,c2||

There can be as many lines as you like for different orderings, however, there has to be only two columns, both of which should contain all names for layers of interest.

Similar to the utilization of *samples information* file, you need to create a samples database for this file as well. To make things simpler, I will use the previous *samples information* file, along with this new order file to create a new samples database prior to calling anvi-interactive:

{% highlight bash %}
anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o samples.db
anvi-interactive -t tree.txt -d view_data.txt -p new_profile.db --title 'Interface Demo VII' --manual -s samples.db
{% endhighlight %}

In this new display you will find out that the two orderings (`test_tree` and `test_list`) appears in the '**Sample order**' combo box under in the '**Samples**' tab on the left panel, along with all the orderings autmoatically generated based on the samples information file:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/07-left.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/07-left.png" /></a>
</div>

Selecting the `test_tree` order, and re-drawing the tree will result in a tiny dendrogram for the layers of interst:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/07.png"><img src="{{ site.url }}/images/anvio/2016-02-27-the-anvio-interactive-interface/07.png" width="50%" /></a>
</div>

{:.notice}
Anvi'server link (don't forget to play with the Samples order combo box): [http://anvi-server.org/public/meren/interface_demo_VII](http://anvi-server.org/public/meren/interface_demo_VII){:target="_blank"}

# 2. Using the anvi'o interactive interface

We tried our best to make anvi'o interface intuitive, easy-to-learn, and easy-to-use. In an ideal world, you shouldn't need a tutorial to start using it, and learn it through trial and error. But here is a very general overview of the major components of the interface.

If you happen to realize there is something missing, or it would have been helpful to you if something was better explained in this document, please do not hesitate to drop us a line.

## 2.1 An overview of the display

The interactive interface has two major areas of interaction: the space for visualization on the right, and the left panel. The left panel gives access to various controls to work with the data visualized, and improve the presentataion of it.

## 2.2 The left panel

At the bottom of the layers tab there is a section with tiny controls that are available in all tabs. Through these controls you can,

* **Create** or **refresh** the display when necessary using the draw button (some changes require you to do that),
* **Zoom in, zoom out, and center** the display.
* **Download your display as an SVG file**. 

{:.notice}
See the post "[Working with SVG files anvi'o generate]({% post_url anvio/2016-10-27-high-resolution-figures %})" for tips about how to work with large SVG files using Inkscape.

### 2.2.1 Layers tab

Through the layers tab you can,

* **Change general settings for the tree** (i.e., switching between circle or rectengular displays, changing tree radius or width), **and layers** (i.e., editing layer margins, or activating custom layer margins).
* **Load or save states** to store all visual settings, or load a previously saved state.
* **Customize individual layers** by switching between different **display modes** depending on the layer type (i.e., 'text' or 'color' mode for categorical layers, or 'bar' or 'intensity' mode for numerical layers), **set normalization** (i.e., 'square-root', or 'log' normalization), **minimum, and maximum** cutoff values for numerical layers, or set **layer height**, and **layer margin** (i.e., its distance from the previous layer).
* Use the **multi-selector** at the bottom to change settings for multiple layers at once.


### 2.2.2 Samples tab

Samples tab is for the additional data you provide the interface through a samples database (see samples order and samples infomration sections above). Through this layer you can,

* **Change the order of layers** using automatically-generated or user-provided orders of layers using the Sample order combo box,
* **Customize individual samples information entries**.

Changes in this tab can be reflected to the current display without re-drawing the entire tree unless the sample order is changed.

### 2.2.3 Bins tab

Anvi'o allows you to create selections of items shown in the display (whether they are metagenomic contigs, 16S rRNA tags, or any other type of information). Bins tab allow you to maintain these selections. Any selection on the tree will be added to active bin in this tab (the state radio button next to a bin defines its activity). Through this tab you can,

* **Create or delete bins**, **set bin names**, **change the color of a given bin**, or sort bins based on their name, the number of units they carry, or completion and contamination estimates (completion / contamination estimates are only computed for genomic or metagenomic analyses).
* View **the number of selected units** in a given bin, and see the **list of names in the selection** by clicking the button that shows the number of units described in the bin.
* **Store a collection of bins**, or **load a previously stored collection**.

### 2.2.4 Mouse tab

The mouse tab displays the value of items underneath the mouse pointer while the user browse the tree.

Displaying the numerical or categorical value of an item shown on the tree is not an easy task. We originally thought that displaying pop-up windows would solve it, but besides the great overhead, it often became a nuisance while browsing parts of the tree. We could show those pop-up displays only when use clicks on the tree, however click-behavior is much more appropriate to add or remove individual items from a bin, hence, it wasn't the best solution either. So we came up with the 'mouse tab'. You have a better idea? I am not surprised! We would love to try improve your experience: please enter an issue, and let's discuss.

### 2.2.5 Search tab

It does what the name suggests. Using this tab you can,

* **Build expressions to search items** visualized in the main display.
* **Highlight matches**, and **append** them to, or **remove** them from the **selected bin** in the Bins tab.  

## Tips and tricks 

Here are some small conveniences that may help the interface serve you better (we are happy to expand these little tricks with your suggestions).

* You can zoom to a section of the display by making a rectangular selection of the area **while the pressing the `shift` button**.

* You can click an entire branch to add items into the selected bin, and remove them by **right-clicking** a branch.

* If you click a branch **while pressing the `Command` or `CTRL` button**, it will create a new bin, and add the content of the selection into that bin.

* By pressing `1`, `2`, `3`, `4`, and `5`, you can go between Layers, Bins, Samples, Mouse, and Search tabs!
