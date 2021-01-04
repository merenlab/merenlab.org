---
layout: page
title: anvi-rename-bins [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Rename all bins in a given collection (so they have pretty names).

See **[program help menu](../../../../vignette#anvi-rename-bins)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program **creates a new <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> from the <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span>s in another collection with specific guidelines.** This is especially helpful when you want to merge multiple collections later or share your project with someone, and you want all of your bins to have nicer names than the default `bin_01`, `bin_02`, etc. based on the order you binned them in. 

So let's take a look at what this program can do with a simple example. 

### Example 1: Renaming all bins in a collection 

Let's say you have a collection called `MY_COLLECTION`, which has four bins: `really`, `bad`, `bin`, and `names`. These names just won't do, so let's get to renaming. To rename all of my bins and put them into a collection called `SURFACE_OCEAN_SAMPLES`, you could run 

<div class="codeblock" markdown="1">
anvi&#45;rename&#45;bins &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                 &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                 &#45;&#45;prefix SURFACE_OCEAN \
                 &#45;&#45;collection&#45;to&#45;read MY_COLLECTION \
                 &#45;&#45;collection&#45;to&#45;write SURFACE_OCEAN_SAMPLES \
                 &#45;&#45;report&#45;file rename.txt
</div>

And voila! Now you have a second collection named `SURFACE_OCEAN_SAMPLES` that contains your four bins, now named  `SURFACE_OCEAN_Bin_00001`, `SURFACE_OCEAN_Bin_00002`, `SURFACE_OCEAN_Bin_00003`, and `SURFACE_OCEAN_Bin_00004`. The order that the numbers are in represents the quality of the bin as a MAG, given by the completion minus redunancy. 

The file `rename.txt` is just a tab-delimited file that contains a summary of your renaming process. The first column has the original name of the bins that you renamed, the second has their new names, and the remaining columns contain information about those bins (like their completion, redundency, and size). 

### Example 2: Separating out the MAGs 

Okay, but what if you want to label your MAGs separately from your bins? You don't like `SURFACE_OCEAN_bin_00004` since it only has a completition stat of 50 percent, and you're not sure if you want to include `SURFACE_OCEAN_bin_00003`  since it has 50 percent redundency. How can you differenciate these iffy bins in your collection? 

Here is the solution: 

<div class="codeblock" markdown="1">
anvi&#45;rename&#45;bins &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                 &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                 &#45;&#45;prefix SURFACE_OCEAN \
                 &#45;&#45;collection&#45;to&#45;read MY_COLLECTION \
                 &#45;&#45;collection&#45;to&#45;write SURFACE_OCEAN_MAGS \
                 &#45;&#45;report&#45;file rename.txt \ 
                 &#45;&#45;call&#45;MAGs \
                 &#45;&#45;min&#45;completition&#45;for&#45;MAG 70 
</div>

Now, the collection `SURFACE_OCEAN_MAGS` will include  `SURFACE_OCEAN_MAG_00001`, `SURFACE_OCEAN_MAG_00002`, `SURFACE_OCEAN_MAG_00003`, and `SURFACE_OCEAN_Bin_00004`. These are exactly the same bins that the collection contained before, but now the names differenciate the wheat from the chaff. 

Now, let's make that same collection (still called `SURFACE_OCEAN_MAGS`) that doesn't include `SURFACE_OCEAN_Bin_00003` as a MAG, since the redundency is too high for what we want to look at right now. 

<div class="codeblock" markdown="1">
anvi&#45;rename&#45;bins &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                 &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                 &#45;&#45;prefix SURFACE_OCEAN \
                 &#45;&#45;collection&#45;to&#45;read MY_COLLECTION \
                 &#45;&#45;collection&#45;to&#45;write SURFACE_OCEAN_MAGS \
                 &#45;&#45;report&#45;file rename.txt \ 
                 &#45;&#45;min&#45;completition&#45;for&#45;MAG 70 \
                 &#45;&#45;max&#45;redundancy&#45;for&#45;MAG 30 \
                 &#45;&#45;call&#45;MAGs
</div>

Now `SURFACE_OCEAN_MAGS`   will include  `SURFACE_OCEAN_MAG_00001`  `SURFACE_OCEAN_MAG_00002`,  `SURFACE_OCEAN_Bin_00003`, and `SURFACE_OCEAN_Bin_00004`.

You also have the option to only classify bins above a certain minimum size as MAGs. 

### Example 3: An example use case in a workflow

For an example use case, on [this page](http://merenlab.org/tutorials/infant-gut/#renaming-bins-in-your-collection-from-chaos-to-order), anvi-rename-bins is used to create a new collection called `MAGs` that contains differenciates bins that have a completion stat of more than 70 percent, and renames all of those bins with the prefix `IGD` (which stands for infant gut dataset). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-rename-bins.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-rename-bins) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
