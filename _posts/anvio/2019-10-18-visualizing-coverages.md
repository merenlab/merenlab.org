---
layout: post
title: "Visualizing contig coverages to make informed statements regarding microbial population diversity"
excerpt: "How to take mapping results, visualize them in anvi'o, and export them for publication"
modified: 2019-10-18
categories: [anvio]
comments: true
authors: [emily, ryan]
---

{% include _toc.html %}

{: .notice}
This workflow is for `v6` and later versions of anvi'o. You can identify which version you have on your computer by typing `anvi-self-test --version` in your terminal.

{% capture images %}{{site.url}}/images/anvio/2019-10-18-visualizing-coverages{% endcapture %}


"Read recruitment", "mapping reads" or even just "mapping" are words that get thrown around all the time. And if you don't know what they mean, this process can be a little confusing. But it can also be extrememly powerful. If you already feel confident with read recruitment, feel free to skip down to the new stuff, probably starting at `What is this tutorial about?`. If you want a refresher, or a slighly more in-depth understanding, here's a quick summary about the process and the utility of read mapping (aka read recruitment, aka mapping - they all mean the same thing).

# How *does* read recruitment work?

Understanding superficially how read mapping works is fairly straighforward. You have a genome(s) or assembled contigs and you're figuring out how many of your short reads from your sequenced samples match to this reference. Like the figure below shows, you can do this for multiple metagenomes (ie. multiple samples) to compare their similarities and differences. It also gives you important information like the coverage (the mean nubmer of reads that map to that location) or the detection (how much of your reference is covered by at least one read).

[![coverages]({{images}}/read_mapping_figure.png)]({{images}}/read_mapping_figure.png){:.center-img .width-90}

Intuitively, you can think of taking each read, determining where and if this sequence matches to the reference, and placing that read there. That's how the data is displayed in coverage plots, and often this level of understanding is good enough for our purposes. But this is not exactly how it works. Although theoretically this approach would work, it is incredibly computationally intensive and would probably take years (depending on the size of the reference and the number of reads) to run most of the mapping experiments that we'd like to do. Instead, read recruitment software uses some fancy approaches to exponentially increase the speed of mapping.

If you're interested, here is a short summary (to the best of our understanding) of how a mapping index is made, and how the reads are aligned. This is taken almost completely from this [video](link) and the original [bowtie paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2690996/). The slides that accompary the video can be found [here](https://ocw.mit.edu/courses/biology/7-91j-foundations-of-computational-and-systems-biology-spring-2014/lecture-slides/MIT7_91JS14_Lecture5.pdf). This is meant as a short introduction to the complexities of mapping, and by no means a full explanation. If you're interested in learning more in depth how this process works, see one of the afore mentioned references.

In our lab we use bowtie2 to recruit reads, so this explanation will be in the context of bowtie, although most read mapping software use similar techniques. Typically, you start with a reference file input in FASTA format. For simiplicity, let's say this is a reference genome (although this could also be contigs from an assembly). Bowtie takes this FASTA file and uses a Burrows-Wheeler Transformation to convert it to a bowtie index. A Burrows-Wheeler transformation is a method of compressing your data by cyclically rotating the sequence to get all possible rotations of the sequence, ordering them alphabetically, ranking each letter by its occurance, and keeping only the right most column of data. The ordering means that the rank, the nth time that that letter has appeared in the sequence, is the same in the first and last columns. The letters are also assigned a count value, in other words their original position in the sequence. We can use the Last to First (LF) function to match a character with its rank (these values are stored by the algorithm) and use the Walk Left function to reconstruct the initial sequence from these data. More information on these functions can be found [here](https://ocw.mit.edu/courses/biology/7-91j-foundations-of-computational-and-systems-biology-spring-2014/lecture-slides/MIT7_91JS14_Lecture5.pdf).

Ok, so now we have a bowtie index, how do we actually map the reads? Our input reads will typically be in FASTQ format, so we'll have the read information and a quality score for each base. To find reads that match exactly, we use something called a Full-text index in Minute space (FM) index. This matches your read from right to left (5' to 3') to the bowtie index. Each iteration of top and bottom using the FM index indicates the range of rows in the index that have progressively longer matches to our query read. If the range is 0, it means that there's nothing that matches the read sequence exactly and it wouldn't be mapped. But we can immediately think of a problem here. What if there are sequencing errors in our data, or more importantly, true biological variation? We want to be able to map reads that aren't a perfect match. To get around this, bowtie uses a backtracking alignment search. When it encounters a mismatch, it attempts to add all of the differernt bases to see what matches and then continues. It can also try to backtrack to the lowest quality match, replace this with another base and continue from there. The number of mismatches and amount of backtracking is limited to avoid reads mapping that do not belong in that position. This algorithm is "greedy", meaning it will keep the first mismatch that works but not necessarily the best one. Bowtie can be forced to use the best one by supplying the flag `--best`. Because bowtie restricts the amount of backtracking, it may seem like it is biased toward the right hand side of the read, however the left hand side of the read is typically higher quality, and bowtie can use a mirrow index to avoid this bias. A mirror index is the reference sequence reverse complementd and made into a bowtie index, then mapping reversed reads to this index. Again, see [here](https://ocw.mit.edu/courses/biology/7-91j-foundations-of-computational-and-systems-biology-spring-2014/lecture-slides/MIT7_91JS14_Lecture5.pdf) for more information.

Alright, the reads have now been mapped, but where are these reads in the reference? To determine this, we take the suffix array, or in other words the file that was generated that kept track of the postition of every base in the original sequence. Unfortunately, this complete file is very large, so instead of keeping the entire thing, bowtie keeps every 32nd entry (this number can be adjusted) and then uses the Walk-Left function (which reconstucts the initial sequence) until it hits the suffix of interest.

Yay! The reads are mapped. Now on to the biologically relevant part of the post.

## Why are we looking at coverage plots?

We often use read recruitment results to make estimates of relative bacterial abundance. For any given section of reference DNA, we can tell what percentage of the reads from that metagenome aligned to this reference. In anvi'o, we can show exactly where the reads mapped and how many mapped to each nucleotide position by visualizing coverage plots from the inspect page (see below). This is often more helpful than just knowing the mean coverage of that reference because we can tell if the entire reference is covered evenly, or if all the reads are being recruited to one section, thus inflating the mean coverage. Although, if most reads are mapping to one section of the contig this can be biologically interesting too!

For example, here is an example of a contig from an assembly of a human gut metagneome:

[![coverages]({{images}}/ITA0013_cov_plot_blog.png)]({{images}}/ITA0013_cov_plot_blog.png){:.center-img .width-90}

 The mean coverage of this contig is probably  around 25, but even a quick glance at this tells us that the number 25 is pretty meaningless. The depth of coverage of the first section of this contig is ~100X, but the rest of the contig looks like it has ~10X coverage at best. If we take a mean coverage value of 25 it suggests that all the genes in this contig are present uniformly, which is clearly not the case.

In addition to coverage, we can also use these results to visualize the diversity within a sample by looking at single nucleotide variants (SNVs). A single nucleotide variant occurs when a nucleotide in the mapped read is different than the reference (these are the mismatches that were discussed in the section about how bowtie works). In the image above, the SNVs are the green, red and black bars extending up from the x-axis.  Green indicates a synonomous mutation, red is non-synonomous, and black is a SNV in a non-coding region. The height of the bar indicates the proportion of reads that have a different nucleotide than the reference. Visualizing SNVs can be very powerful if you want to make statements about monoclonality of a sample, or to identify potential nucleotides that are under selective pressures in that environment (see this [post](http://merenlab.org/2018/09/04/structural-biology-with-anvio/) for more information).

*If you want to know more about mapping, check out the genome-resolved metagenomics section of the [infant gut](http://merenlab.org/tutorials/infant-gut/#chapter-i-genome-resolved-metagenomics) tutorial or the the blog post on [coverage variation](http://merenlab.org/2016/12/14/coverage-variation/).*

## What is this tutorial about?

The purpose of this tutorial is not to convince you that manual inspection of your contigs is well worth your time (although maybe now you're a little more convinced), but rather to explain new features in anvi'o v6 that make this process managable - **especially** if you have a ton of samples.

If you have any questions, please get in touch with us and/or other anvians:

{% include _join-anvio-slack.html %}

This post assumes that you have profile, auxiliary-data and contigs databases for one or more samples mapped to one or more contigs. For a more in depth explanation of this process please look [here](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). The post you're currently on will walk you through the process of inspecting your splits using [anvi-inspect](http://merenlab.org/software/anvio/vignette/#anvi-inspect) for a quick inspection and [anvi-script-visualize-split-coverages](insertlink) to export a PDF of these results. To make things more interesting, we'll do a re-analysis of publically available data mapped to a globally prevalent marine plasmid from a recent [study](https://www.pnas.org/content/early/2019/09/19/1905878116). If you'd like to download the `PROFILE.db` `AUXILIARY-DATA.db` and `CONTIGS.db` to follow along with this tutorial, go [here](add_in_link).

{:.notice}
A quick reminder definition: anvi'o will divide contigs longer than 25,000 basepairs into multiple "splits" to visualize them in the interactive interface. If you only have one contig, as is our case, it can still be referred to as a single split.

{:.notice}
It's not strictly necessary for this tutorial, but the brief rationale for using these particular data is as follows. From marine metagenomes, [Petersen et al](https://www.pnas.org/content/early/2019/09/19/1905878116) identified a highly conserved and globally distributed plasmid called pLA6\_12. The backbone of the plasmid is conserved, but there is a hyper variable section that can have variable genes, such as chromate resistance, that appear to benefit the bacterial cell. The authors claim that the backbone of this plasmid backbone is conserved between sampling locations but don't explain how they came to this conclusion, nor do they give any details in the methods section of the paper. We can use read recruitment of metagenomes to this plasmid to examine their claims. We downloaded some of the metagenomes that they specify in the supplemental data, mapped these to pLA6\_12, and built a `PROFILE.db` and `CONTIGS.db`. In this tutorial we'll visualize these.

## How to visualize coverage plots directly from the interactive interface

Before jumping into `anvi-inspect` and `anvi-script-visualize-split-coverages` let's go through the tried and true approach of visualizing coverage plots once you've opened the anvi'o interface. This will give you a better idea of the utility of `anvi-inspect` and `anvi-script-visualize-split-coverages`. If you're comfortable with inspecting contigs from the interface, feel free to skip this section.

Take your `PROFILE.db` and `CONTIGS.db` and open the interactive interface:

```bash
anvi-interactive -p PROFILE.db -c CONTIGS.db
```

This produces a figure that looks like this. It's certainly not the prettiest anvi'o figure you've ever seen, but we just care about the inspect page, so let's move on.

[![coverages]({{images}}/blog_interactive_interface.png)]({{images}}/blog_interactive_interface.png){:.center-img .width-90}

We only have one contig (pLA6_12), so we can right-click anywhere and click inspect to open the inspect page that looks like this (below). If you had more than one contig, you'd just pick the contig you want to inspect and do the same thing. This is great, but sometimes you'll have a lot of samples, and opening the interactive interface and then the inspect page will take a **LONG** time.

[![coverages]({{images}}/ScreenShot_inspect_pLA6_12.png)]({{images}}/ScreenShot_inspect_pLA6_12.png){:.center-img .width-90}


## How to use `anvi-inspect` to visualize coverage plots directly

To avoid having to open the interactive interface and instead going straight to the inspect page, you can use the command `anvi-inspect`. For this you'll need the same `PROFILE.db` and `CONTIGS.db`, but you'll have to specify the split that you're interested in. For example:

```bash
anvi-inspect -p PROFILE.db \
             -c CONTIGS.db \
             --split-name pLA6_12_000000000001_split_00001
```

will give you the same inspect page that we saw before, but allows you to skip going through the interactive interface (if you don't know what your split name is, just keep reading and we'll show you in a second how to find that out).

Ok, great. We have some nice looking split coverages. Let's say that now we wanted to export these coverage plots. We could take screen shots of the inspect page but that's not very elegant, nor is it reproducible. And it's even less ideal if you have a lot of samples.

## Exporting coverage plots to a PDF

Thanks to Ryan Moore at the University of Delaware (**insert his github and twitter here**), we now have a program that will take inputs from one or more files to produce these coverage plots as a PDF. In other words, this program will give you the inspect page produced in `ggplot2`. If you don't already have it, you'll need to download R and install the `ggplot2` and `optparse` packages.

### I just want the coverage of my split across samples:

The most basic funcitonality of `anvi-visualize-split-coverages` is to export coverage plots as a PDF. In order to do this, it needs to know the coverage of your split. To get this, you can run the command [anvi-get-split-coverages](http://merenlab.org/software/anvio/vignette/#anvi-get-split-coverages) and specify the split name similar to what's shown below.

If you don't know what your split of interest is called, you can run this first:

```bash
anvi-get-split-coverages -p PROFILE.db \
                         -c CONTIGS.db \
                         --list-splits
```

Then to get the split coverages run this:

```bash
anvi-get-split-coverages -p PROFILE.db -c CONTIGS.db \
                         --split-name pLA6_12_000000000001_split_00001 -o pLA6_12_000000000001_split_00001_coverage.txt
```


Now you use the `pLA6_12_000000000001_split_00001_coverage.txt ` output in the `anvi-visualize-split-coverages` command like this:

```bash
anvi-script-visualize-split-coverages -i pLA6_12_000000000001_split_00001_coverage.txt -o pLA6_12_000000000001_split_00001_inspect.pdf
```

[![coverages]({{images}}/default_cov_plot.png)]({{images}}/default_cov_plot.png){:.center-img .width-90}

### I want to generate coverage plots with single nuleotide variants (SNVs):

`anvi-visualize-split-coverages` can plot the coverage values with the corresponding SNVs if you also provide the output of [anvi-gen-variability-profile](http://merenlab.org/software/anvio/vignette/#anvi-gen-variability-profile). SNVs are nice to have if you want to make statements about sub-populations or identify nucletides that are under selective pressures. If you want to learn more about `anvi-gen-variability-profile`, have a look at this [post](http://merenlab.org/2015/07/20/analyzing-variability/#the-output-matrix). If you just want to get to the SNVs on your coverage plots, go straight to running a command like this:

```bash
anvi-gen-variability-profile -p PROFILE.db -c CONTIGS.db \
                             --splits-of-interest split_of_interest.txt \
                             --include-split-names \
                             -o pLA6_12_000000000001_split_00001_SNVs.txt \
                             --include-contig-names
```

Make sure you include the `--include-contig-names` and `--include-split-names` flags. The `split_of_interest.txt` file is the split you're interested in a plain text file, like this:

```
pLA6_12_000000000001_split_00001
```

Then run the `anvi-visualize-split-coverages` like this:

```bash
anvi-script-visualize-split-coverages -i pLA6_12_000000000001_split_00001_coverage.txt \
                                      -o pLA6_12_000000000001_split_00001_inspect.pdf \
                                      --snv-data pLA6_12_000000000001_split_00001_SNVs.txt
```

[![coverages]({{images}}/cov_plot_SNVs_pla6_12.png)]({{images}}/cov_plot_SNVs_pla6_12.png){:.center-img .width-90}

Yay! Now we have SNVs. The default SNV colors are green, red and grey.

* Green = synonomous (the amino acid doesn't change)
* Red = non-synonomous (the amino acid does change)
* Grey = intergenic region SNV.

As we can see from the plot above, the coverage in the above samples varies widely, and while this is good to distinguish highly covered vs lowly covered samples, it also makes it difficult to visualize the ones with less coverage, especially of you want to look at SNVs. Thankfully, there is a solution for this. You can run exactly the same `anvi-visualize-split-coverages` command, but add in the flag `--free-y-scale TRUE`.

```bash
anvi-script-visualize-split-coverages -i pLA6_12_000000000001_split_00001_coverage.txt \
                                      -o pLA6_12_000000000001_split_00001_inspect.pdf \
                                      --snv-data pLA6_12_000000000001_split_00001_SNVs.txt
                                      --free-y-scale TRUE
```

This will give you a plot that looks like this:

[![coverages]({{images}}/free_y_scale.png)]({{images}}/free_y_scale.png){:.center-img .width-90}

### I have a lot of samples and I don't want to crash ggplot

Sometimes you'll have so many samples that it may crash ggplot trying to plot them all together. Or maybe you want separate PDFs for each of the sets of samples. Or maybe you just like messing around with new anvi'o commands. Regardless of your motivation, there is a way to do this.

Another additional file you can create is a `samples_data.txt` file. This will specify which samples should be grouped together into the same PDF. For example, some of the metagenomes we downloaded were from the Red Sea, and some were metagenomes from pelagic zones all over the world (the Malaspina metagenomes). If I wanted to Red Sea mapping results to one PDF, and the Malaspina metagenomes mapping to another PDF, I would create a file that looks like this:

```
sample_name	sample_group
MSP0109	MAL
MSP0112	MAL
MSP0114	MAL
MSP0121	MAL
MSP0131	MAL
MSP0144	MAL
MSP0146	MAL
RED0041	RED
RED0045	RED
```

Then I would run `anvi-visualize-split-coverages` like this:

```bash
anvi-script-visualize-split-coverages -i pLA6_12_000000000001_split_00001_coverage.txt \
                                      -o pLA6_12_000000000001_split_00001_inspect.pdf \
                                      --snv-data pLA6_12_000000000001_split_00001_SNVs.txt \
                                      --sample-data sample_data.txt
```

And it will give me two separate PDFs, each with those samples.

#### Malaspina:

[![coverages]({{images}}/MAL_plots.png)]({{images}}/MAL_plots.png){:.center-img .width-90}

#### Red sea:

[![coverages]({{images}}/RED_plots.png)]({{images}}/RED_plots.png){:.center-img .width-90}


### I want prettier and customizable coverage plots

The outputs above are great, but they just give you black coverage plots and default SNV colors. Plus, everything has the same y-axis even if the coverage values are radically different. But of course these can be changed to fit the users needs. For example, let's say I wanted the sample plots we outputed above, but I wanted to color them based on where the samples are from, restrict the maximum coverage to 2000, allow a sliding y-axes for each plot, and have customized SNV lines I would run:

```bash
anvi-script-visualize-split-coverages -i pLA6_12_000000000001_split_00001_coverage.txt \
                                      -o pLA6_12_000000000001_split_00001_inspect.pdf \
                                      --sample-data sample_data_colors.txt \
                                      --snv-data pLA6_12_000000000001_split_00001_SNVs.txt \
                                      --free-y-scale TRUE \
                                      --max-coverage 2000 \
                                      --snv-marker-transparency 0.9 \
                                      --snv-marker-width 0.1
```

Where the `sample_data_colors.txt` file is:

```
sample_name	sample_group
MSP0109	MAL	#000080
MSP0112	MAL	#000080
MSP0114	MAL	#000080
MSP0121	MAL	#000080
MSP0131	MAL	#000080
MSP0144	MAL	#000080
MSP0146	MAL	#000080
RED0041	RED	#8b0000
RED0045	RED	#8b0000
```

And voila!

#### Malaspina:

[![coverages]({{images}}/MAL_color_plots.png)]({{images}}/MAL_color_plots.png){:.center-img .width-90}

#### Red sea:

[![coverages]({{images}}/RED_color_plots.png)]({{images}}/RED_color_plots.png){:.center-img .width-90}

This is of course just one example of a customizable plot and the full list of optional variations is available in the help menu for `anvi-script-visualize-split-coverages` (`anvi-script-visualize-split-coverages -h`).

Now we've visualized the mapping results from [Petersen et al](https://www.pnas.org/content/early/2019/09/19/1905878116), we can try to make some fair conclusions. In their paper, they claim that this plasmid is nearly identical between all of these metagenomes, aside from the variable region that we can see in the middle. Our mapping results confirm that this is indeed the case. There are a lot of SNVs in comparison to the reference pLA6_12 plasmid, but clearly it is conserved across these environments. Perfect. They could have given more details about how they came to these conclusions, but their conclusions seem valid.
However, if they had looked at their coverage plots, they would have realized that in samples where there appars to be two versions of the plasmid (eg. MAP0144 and MSPO0146 where some reads are mapping to the hypervariable region) the SNV bars do not extend completely to the top of the plot. This indicates that the version that has the reads mapping to the hypervariable region may also have a slightly different backbone that is more similar to the reference than the rest of the environment. But at this point it's all just speculation from us.

Thank you very much if you made it all the way to the bottom and if you still have questions about `anvi-inspect` or `anvi-script-visualize-split-coverages` feel free to reach out to us.

{% include _join-anvio-slack.html %}

Happy inspecting!



