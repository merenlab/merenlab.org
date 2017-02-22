---
layout: page
title: A tutorial on assembly-based metagenomics
modified: 2016-07-13
excerpt: "By Meren Lab and anvi'o developers"
comments: true
image:
  feature: feature-infant-gut-tutorial.png
---

{% include _toc.html %}

{:.notice}
The version number of this tutorial is `1.0`, and for now it is tailored for Illumina paired-end shotgun sequencing with large inserts (i.e., no substantial overlap between two reads in a given pair).

{:.notice}
This is our **very initial attempt** to put together a comprehensive tutorial. If you would like to change something on this page, you can directly **[edit its source code](https://github.com/meren/meren.github.io/blob/master/tutorials/assembly-based-metagenomics/index.md)** by clicking the _"Edit this file"_ icon on the right top corner of the page, and send us a pull request through GitHub. We will be very thankful. You will need to login to GitHub to see that icon.

**Assembly** and **mapping** are key steps for most **assembly-based, genome-resolved metagenomic studies**, and there are many ways to accomplish each of these steps. That's why, the [anvi'o metagenomic workflow]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}) only starts once you have your contigs and BAM files available.

But a metagenomic study starts much earlier than assembly and mapping. Experimental design, sampling and storage strategies, library preparation, and sequencing are all critical steps that should be given considerable thought, especially if the purpose is to recover metagenome-assembled genomes.

Despite the diversity of ways to get your contigs and BAM files, or preparing your libraries and perform your sequencing, we thought we could create a resource to explain the metagenomic workflow [we]({{site.url}}/people/) and our colleagues have been using in our lab and at the [MBL](http://mbl.edu).

This tutorial will take you from notes on sampling and library preparation considerations for sequencing, to assembled contigs and BAM files, at which point you will be ready to follow the anvi'o metagenomic workflow, or any other platform to make sense of your data.


# Before Sequencing

{:.notice}
You have the question.

This section mostly deals with the field and wet-lab aspect of the workflow.

## Sampling

Storage? Replicates?

## Library preparation

DNA extraction? Size selection?

## Sequencing

Which sequencing platform? How many samples per lane?

# After Sequencing

{:.notice}
You have your sequences.

This section mostly deals with the computational aspect of the workflow.

If you have some familiarity with the UNIX shell environment, you will have no trouble following the tutorial with your own data. We are working on a more automated workflow which should help people who have no familiarity with the UNIX shell environment, but we need a bit more time for that.

To make things very simple, we will assume you have two samples throughout this tutorial, `Sample_01` and `Sample_02`.

Just to make sure we are mostly on the same page and you have all the software you need installed on your system, here are a couple of commands, and their outputs:

{% highlight bash %}
meren SSH://MBL ~ $ bowtie2 --version | head -n 1 | awk '{print $3}'
2.2.9
meren SSH://MBL ~ $ iu-filter-quality-minoche -v
Illumina-utils v1.4.6
meren SSH://MBL ~ $ samtools --version
samtools 1.3
Using htslib 1.3
meren SSH://MBL ~ $ megahit -v
MEGAHIT v1.0.6
meren SSH://MBL ~ $ anvi-init-bam -v
Anvio version ................................: 2.0.1
Profile DB version ...........................: 16
Contigs DB version ...........................: 6
Samples information DB version ...............: 2
Auxiliary HDF5 DB version ....................: 1
Users DB version (for anvi-server) ...........: 1
{% endhighlight %}

You may be lacking some of these software. You can install anvi'o using methods described in [this article]({% post_url anvio/2016-06-26-installation-v2 %}), or you can take a look at the article that gives [recipes to install third-party software]({% post_url anvio/2016-06-18-installing-third-party-software %}) (what you are looking for may not be covered in that article if it is very simple to install, in which case Google will help you find out how to install those on your system).

Let's start.


## Demultiplexing

{:.notice}
You have `R1`, `R2`, and `I1` files for an entire Illumina HiSeq or MiSeq lane, and you need to demultiplex your samples using some barcodes.

The TAB-delimited `barcodes_to_samples.txt` file should look like this:

|Sample_01|TAATCGGATTCC|
|Sample_02|GCACACACGTTA|

where each line has two columns to describe the sample name and its barcode.

If everything is ready, these commands will demultiplex your samples:

{% highlight bash %}
meren SSH://MBL ~ $ mkdir 00_RAW/
meren SSH://MBL ~ $ iu-demultiplex -s barcode_to_sample.txt --r1 R1.fastq --r2 R2.fastq --index I1.fastq -o 00_RAW/
meren SSH://MBL ~ $ ls 00_RAW/
00_DEMULTIPLEXING_REPORT
Sample_01-R1.fastq     
Sample_01-R2.fastq
Sample_02-R1.fastq 
Sample_02-R2.fastq
{% endhighlight %}

It is always a good idea to take a quick look at the report file to make sure everything seems alright:

{% highlight bash %}
meren SSH://MBL ~ $ column -t 00_RAW/00_DEMULTIPLEXING_REPORT
{% endhighlight %}

Now you have everything you need and you can continue with quality filtering.

## Quality Filtering

{:.notice}
You have raw `R1` and `R2` files for `Sample_01` and `Sample_02`, and you need to do quality filtering.

You first need to generate a TAB-delimited `samples.txt` file to point out where are your raw `R1` and `R2` files for each sample:

|sample|r1|r2|
|:--|:--:|:--:|
|Sample_01|00_RAW/Sample_01-R1.fastq|00_RAW/Sample_01-R2.fastq|
|Sample_02|00_RAW/Sample_02-R1.fastq|00_RAW/Sample_02-R2.fastq|


Then you need create a directory for quality-filtered `R1` and `R2`, and then use `iu-gen-configs` program with `samples.txt` to crate config files for illumina-utils in it:

{% highlight bash %}
meren SSH://MBL ~ $ mkdir 01_QC
meren SSH://MBL ~ $ iu-gen-configs samples.txt -o 01_QC
meren SSH://MBL ~ $ ls 01_QC/
Sample_01.ini
Sample_02.ini
{% endhighlight %}

Now you are ready to run quality filtering for each of your samples. You can do it for one of them the following way:

{% highlight bash %}
meren SSH://MBL ~ $ iu-filter-quality-minoche 01_QC/Sample_01.ini
{% endhighlight %}

{:.notice}
You should use `iu-filter-quality-minoche` only if you have large inserts. If you have partially overlapping reads, you should use `iu-merge-pairs` program at this step.

Alternatively you can do it for all your samples at once in a `for` loop, instead of doing it one by one:

{% highlight bash %}
meren SSH://MBL ~ $ for ini in 01_QC/*.ini; do iu-filter-quality-minoche 01_QC/$ini; done
{% endhighlight %}

{:.notice}
Of course, if you have access to a cluster, you should add necessary modifications to these commands.

Once you are done, the contents of the `01_QC/` directory should look like this:

{% highlight bash %}
meren SSH://MBL ~ $ ls 01_QC/
Sample_01-QUALITY_PASSED_R1.fastq
Sample_01-QUALITY_PASSED_R2.fastq
Sample_01-READ_IDs.cPickle.z
Sample_01-STATS.txt
Sample_01.ini
Sample_02-QUALITY_PASSED_R1.fastq
Sample_02-QUALITY_PASSED_R2.fastq
Sample_02-READ_IDs.cPickle.z
Sample_02-STATS.txt
{% endhighlight %}

You should definitely take a quick look at `*_STATS.txt` files to see whether things went alright. A STATS file will look like this, and you can find more information about it here:

{% highlight bash %}
meren SSH://MBL ~ $ cat 01_QC/Sample_01-STATS.txt
number of pairs analyzed      : 122929
total pairs passed            : 109041 (%88.70 of all pairs)
  total pair_1 trimmed        : 6476 (%5.94 of all passed pairs)
  total pair_2 trimmed        : 9059 (%8.31 of all passed pairs)
total pairs failed            : 13888 (%11.30 of all pairs)
  pairs failed due to pair_1  : 815 (%5.87 of all failed pairs)
  pairs failed due to pair_2  : 12193 (%87.80 of all failed pairs)
  pairs failed due to both    : 880 (%6.34 of all failed pairs)
  FAILED_REASON_P             : 12223 (%88.01 of all failed pairs)
  FAILED_REASON_N             : 38 (%0.27 of all failed pairs)
  FAILED_REASON_C33           : 1627 (%11.72 of all failed pairs)
{% endhighlight %}

You can take a quick look at all samples by running commands like this:

{% highlight bash %}
meren SSH://MBL ~ $ grep 'total pairs passed' 01_QC/*STATS.txt
{% endhighlight %}

{:.notice}
You can find more information about the options, and ways to visualize the quality filtering results [here](https://github.com/meren/illumina-utils/blob/master/README.md).

If all looks good, you are ready for a co-assembly.

## Co-assembly

{:.notice}
You have quality-filtered `R1` and `R2` files for `Sample_01` and `Sample_02`, and you want to get your `contigs.fa`.

{:.notice}
Here we will use MEGAHIT for the assembly, if your lab has a different favorite, please help us expand this tutorial by adding yours. If you want to do it but don't know how, send us an e-mail!

Let's assume quality-filtered FASTA files for `Sample_01` and `Sample_02` are in a directory called `01_QC/`, which looks like this:

{% highlight bash %}
meren SSH://MBL ~ $ ls 01_QC/*fastq
Sample_01-QUALITY_PASSED_R1.fastq
Sample_01-QUALITY_PASSED_R2.fastq
Sample_02-QUALITY_PASSED_R1.fastq
Sample_02-QUALITY_PASSED_R2.fastq
{% endhighlight %}

To make things simpler and to minimize human error, let's create two environment variables:

{% highlight bash %}
meren SSH://MBL ~ $ R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
meren SSH://MBL ~ $ R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
{% endhighlight %}

They will look like this:

{% highlight bash %}
meren SSH://MBL ~ $ echo $R1s
01_QC/Sample_01-QUALITY_PASSED_R1.fastq,01_QC/Sample_02-QUALITY_PASSED_R1.fastq
meren SSH://MBL ~ $ echo $R2s
01_QC/Sample_01-QUALITY_PASSED_R2.fastq,01_QC/Sample_02-QUALITY_PASSED_R2.fastq
{% endhighlight %}

We are ready to run MEGAHIT.

{% highlight bash %}
meren SSH://MBL ~ $ megahit -1 $R1s -2 $R2s --min-contig-len $MIN_CONTIG_SIZE -m 0.85 -o 02_ASSEMBLY/ -t $NUM_THREADS
{% endhighlight %}

You should replace `$MIN_CONTIG_SIZE` and `$NUM_THREADS` with actual values you want. We usually use 1000 for `$MIN_CONTIG_SIZE`, and 40 for `$NUM_THREADS` (which of course depends on the number of CPUs available on a computer system).

Once the co-assembly is done, we first check the log file to take a look at some of the simple stats such as the size of the largest contigs, average contig length, N50, and to make sure things didn't go south during the assembly:

{% highlight bash %}
meren SSH://MBL ~ $ tail 02_ASSEMBLY/log
    [assembler.cpp             : 188]     unitig graph size: 8343857, time for building: 359.859433
    [assembler.cpp             : 206]     Number of bubbles/complex bubbles removed: 854393/135891, Time elapsed(sec): 11.599862
    [assembler.cpp             : 230]     Unitigs removed in excessive pruning: 218649, time: 2.998776
    [assembler.cpp             : 279]     Number of local low depth unitigs removed: 31409, complex bubbles removed: 44857, time: 142.707443
    [assembler.cpp             : 132]     Total length: 3189757381, N50: 1335, Mean: 547, number of contigs: 5824928
    [assembler.cpp             : 133]     Maximum length: 272005
    [utils.h                   : 126]     Real: 1237.1572	user: 65588.0131	sys: 40.3599	maxrss: 18847612
--- [Tue Jun 21 09:50:48 2016] Merging to output final contigs ---
--- [STAT] 787465 contigs, total 1854094345 bp, min 1000 bp, max 272005 bp, avg 2355 bp, N50 2640 bp
--- [Tue Jun 21 09:51:19 2016] ALL DONE. Time elapsed: 36922.054002 seconds ---
{% endhighlight %}

Then we finalize our `contigs.fa` while simplifying contig names in it, and eliminating some of the short contigs if at the same time if necessary:

{% highlight bash %}
meren SSH://MBL ~ $ mkdir 03_CONTIGS
meren SSH://MBL ~ $ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt
{% endhighlight %}

Once this is done, we have our `contigs.fa` under the directory `03_CONTIGS/`. Time to map things!

## Mapping

{:.notice}
You have the quality-filtered `R1` and `R2` files for `Sample_01` and `Sample_02`, as well as your `contigs.fa`, and you want to get BAM files for your samples.

{:.notice}
Here we will use Bowtie2 for mapping, if your lab has a different favorite, please help us expand this tutorial by adding yours.

Let's first create a new directory for mapping results:

{% highlight bash %}
meren SSH://MBL ~ $ mkdir 04_MAPPING
{% endhighlight %}

Next, we need to build an index for our contigs:

{% highlight bash %}
meren SSH://MBL ~ $ bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
{% endhighlight %}

When this is done, this is how you will get an indexed BAM file for one of your samples:

{% highlight bash %}
meren SSH://MBL ~ $ bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 01_QC/Sample_01-QUALITY_PASSED_R1.fastq -2 01_QC/Sample_01-QUALITY_PASSED_R2.fastq -S 04_MAPPING/Sample_01.sam
meren SSH://MBL ~ $ samtools view -F 4 -bS 04_MAPPING/Sample_01.sam > 04_MAPPING/Sample_01-RAW.bam
meren SSH://MBL ~ $ anvi-init-bam 04_MAPPING/Sample_01-RAW.bam -o 04_MAPPING/Sample_01.bam
meren SSH://MBL ~ $ rm 04_MAPPING/Sample_01.sam 04_MAPPING/Sample_01-RAW.bam
{% endhighlight %}

{:.notice}
You should replace `$NUM_THREADS` with an actual integer values you prefer. We usually use 8 for `$NUM_THREADS` when we run a lot of things in parallel (which of course depends on the number of CPUs available on a computer system).

Clearly you will need to do for all your samples. You can take a look at this snippet to see how you can utilize your `samples.txt` file (see the section *Quality Filtering* if you do not have one):

{% gist meren/5c8616959f3eea1f0632b50e2f02fb1e %}

Did it work? Well, we certainly hope it did. Do not forget to make sure you have all your BAM files in the mapping directory, and they have reasonable sizes:

{% highlight bash %}
meren SSH://MBL ~ $ ls -lh 04_MAPPING/
-rw-rw-r-- 1 meren merenlab 476M Jul 14 10:25 Sample_01.bam
-rw-rw-r-- 1 meren merenlab 1.9M Jul 14 10:26 Sample_01.bam.bai
-rw-rw-r-- 1 meren merenlab 772M Jul 14 10:29 Sample_02.bam
-rw-rw-r-- 1 meren merenlab 1.9M Jul 14 10:30 Sample_02.bam.bai
{% endhighlight %}

{:.notice}
At this point the size of each BAM file will be proportional to the number of reads from your metagenomes mapped to your contigs, which is an indication of how well your contigs represent the metagenome.

If you are here, it means you have your `contigs.fa`, and your BAM files! Congratulations! :)

## What is next?

Now you can continue with the [anvi'o metagenomic workflow]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}), or you can take a look at all the other anvi'o stuff using this surprisingly fancy menu:

{% include _project-anvio-navigation.html %}

<div style="margin: 50px;"></div>

Please let us know if something doesn't work.

<div style="margin: 50px;"></div>
