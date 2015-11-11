---
layout: post
title: "Anvi'o Samples Database"
excerpt: "A nifty feature that makes your trees triple awesome."
modified: 2015-11-10
tags: []
categories: [anvio]
comments: true
author: meren
---

{% include _toc.html %}

The anvi'o interactive interface is quite flexible, and any tree can be expanded with additional layers of information. However, until the anvi'o [version 1.1.0](https://github.com/meren/anvio/releases/tag/v1.1.0), it hasn't been so straightforward to add data layers that are relevant to sample layers instead of contigs.

To bridge this gap we developed a new feature that increases the efficacy of the anvi'o's visualization strategy even furhter.

This article will cover everything about this optional feature, **the samples database**.

## The samples database: with and wihtout

Here are two screenshots from the infamous [`mini_test`]({% post_url anvio/2015-05-01-installation %}#running-the-mini-test) just to give a very quick idea about what does the samples database add to an anvi'o analysis.

When you invoke the interactive interface on `mini_test` without a samples database, this is what you see:

{% highlight bash %}
cd anvio/tests/sandbox/test-output
anvi-interactive -p 204-MERGED/PROFILE.db -c CONTIGS.db
{% endhighlight %}

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-10-samples-db/without-samples-db.png"><img src="{{ site.url }}/images/anvio/2015-11-10-samples-db/without-samples-db.png" width="45%" /></a>
</div>

In contrast, when you run the interactive interface on `mini_test` with a samples database, this is what you see (note the new paramter `-s`):

{% highlight bash %}
anvi-interactive -p 204-MERGED/PROFILE.db -c CONTIGS.db -s SAMPLES.db
{% endhighlight %}

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-10-samples-db/with-samples-db.png"><img src="{{ site.url }}/images/anvio/2015-11-10-samples-db/with-samples-db.png" width="70%" /></a>
</div>

The question is, where is this SAMPLES.db coming from.

## Generating a SAMPLES.db

The program `anvi-gen-samples-info-database` helps user to generate a samples database.

There are two essential files that can go into this process:

* The **samples information** file,
* The **samples order** file.

They are both TAB-delimited text files, and as their names suggest, they have slightly different purposes. The samples database can be generated using either of these files, or both of them together.


### The 'samples information' file format

The samples information file contains numerical, categorical, or stacked-bar data for each sample. Here is [the matrix file](https://github.com/meren/anvio/blob/master/tests/sandbox/samples-information.txt) for samples information that is used to generate the samples database that resulted in the second figure in this page looked like this:

|samples|num_reads|percent_mapped_reads|visits|bar1;bar2;bar3|
|---------------|:-------------:|:------:|:--------:|:---------------:|
|s204_6M|1000|95|even|3;3;1|
|s204_7M|5000|93|odd|1;3;3|
|s204_9M|10000|90|odd|3;1;3|

If you are not sure how sample names are seen by anvi'o, you can run this command on your terminal to get this information from your profile database:

{% highlight bash %}
sqlite3 test-output/204-MERGED/PROFILE.db 'select value from self where key = "samples";'
{% endhighlight %}

Easy peasy.

### The 'samples order' file format

The samples order file contains information about different 'orderings' of samples, so the user can access to this information from the interface to arrange all view layers. An order can be a comma separated list of sample names, or it can be a newick tree for the organization of samples. Here is [an example matrix file](https://github.com/meren/anvio/blob/master/tests/sandbox/samples-order.txt) for samples order that exemplifies both newick and list ordering of samples:

|attributes|basic|newick|
|----------|:-------------|:------|
|test_order_01||(s204_6M:0.139971,(s204_7M:9.6205e-10,s204_9M:9.6205e-10)Int3:0.139971);|
|test_order_02|s204_9M,s204_7M,s204_6M||

As you can see, for each *attribute* you can either define a basic, or newick attribute, but not both.

### The command line to generate the database

Once you have either of these files, or both of them ready for your analysis, you can run this command to finally generate the samples database:

{% highlight bash %}
anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o SAMPLES.db
{% endhighlight %}

Anvi'o will tell you if it doesn't like your file(s) for some reason.

{:.notice}
<b>IMPORTANT NOTE for MAC users</b>: When you export a TAB-delimited text file from EXCEL, it will not work properly with the pipeline, since EXCEL fails to insert proper "newline" characters into these files. There are many ways to fix the file, but if you are not familiar with this issue, I would suggest you to download this <a href="http://sourceforge.net/projects/dos2unix/">little command line tool</a> to correct EXCEL's failure to perform one of the simplest tasks imaginable.

## Samples tab

When you run the interface with a samples database, depending on the input file(s) you used to generate the database, new orderings and data layers will appear on the 'samples' tab:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-10-samples-db/samples-tab.png"><img src="{{ site.url }}/images/anvio/2015-11-10-samples-db/samples-tab.png" /></a>
</div>

Because anvi'o cares about its users very much, it will attempt to auto-generate *some* orderings depending on the attributes used in the samples information file. These orders will show up in the 'Samples order' combo box, along with the special orderings the user defined through the samples order file (please take a look at the input files again and orient yourself with the content of the combo box):

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-10-samples-db/sample-order-combo.png"><img src="{{ site.url }}/images/anvio/2015-11-10-samples-db/sample-order-combo.png" /></a>
</div>

{:.notice}
When you change to the order of sample information layers, play with their colors, normalization values, or min/max values, you can simply click the 'Redraw sample layers' button instead of re-drawing the entire tree from scratch. Changing the order, on the other hand, will require the entire tree to be re-drawn.

For instance, if you select a tree-based order from the combo box, and click 'Draw', not only sample layers are going to be ordered according to the structure of the tree, but also a tiny tree will appear on the side:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-10-samples-db/tree-order.png"><img src="{{ site.url }}/images/anvio/2015-11-10-samples-db/tree-order.png" width="60%" /></a>
</div>

Please do not hesitate to ask if you run into any issues with the samples database!

You can visit [http://anvio.org/demo](http://anvio.org/demo) to play with the anvi'o interactive interface running the `mini_test` with a samples database generated using the two files mentioned in this post.