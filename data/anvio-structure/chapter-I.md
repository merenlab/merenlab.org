---
layout: page
title: Chapter I - Reproducing Kiefl et al, 2022
modified: 2021-10-21
excerpt: "A complete reproducible workflow of the manuscript 'FIXME' by Kiefl et al"
comments: true
authors: [evan]
redirect_from:
  - 
---

{% capture images %}{{site.url}}/data/anvio-structure/images{% endcapture %}
{% capture command_style %}background: #D7484822; border: 4px solid #D74848;{% endcapture %}

## Introduction

This multi-part document details the story of how our paper came to be. It includes step-by-step procedures for reproducing our study, in-depth explanations of our analysis decisions, and technical details that were not suitable for the main text or supplemental info.

### Philosophy

A central tenant of the scientific process is reproducibility, yet we are in the midst of a [reproducibility crisis](https://www.nature.com/articles/533452a) in which scientists constantly fail to reproduce published results.

Despite the computer being governed by deterministic rules, computational sciences are not immune to issues of irreproducibility. Only with great care and foresight can a computationally-driven study be considered reproducible. And unfortunately, it is rare to find studies that provide reproducible documents comprehensive enough to replicate even the central findings of the study. Of course this isn't due to ill intent (in the assumed vast majority of cases), but instead due to a lack of time to create such a document, and/or a lack of skills required to make it right.

Indeed, writing this document took a ridiculous amount of time, and structuring my data analysis with reproducibility in mind took several iterations of planning and execution. I have done this analysis not once, but countless times--on different machines and operating systems. And despite all of my efforts, there is tons of room for improvement. This is because the reproducible ceiling is so high. Let's talk more about that.

To some, the holy grail of reproducible workflows is a program script that I'll call `generate_paper.sh`, that when run, downloads all raw data, runs all analyses, creates all figures, tables and numbers, plops them into a LaTeX file where the paper's text is sitting, and compiles a pdf of the paper.

What value does this bring? Well, if a reviewer wants a statistical analysis different from the one used, just open up `generate_paper.sh`, find the relevant section, change `pearson` to `spearman`, and regenerate the paper by typing `./generate_paper.sh`. Or if a reader is intrigued by the implications of a result, they can run `./generate_paper.sh`, find the intermediate files, and extend the analysis to their liking.

The value of `generate_paper.sh` is therefore understandably enormous. But there's a missing element to this, isn't there? It's the _human_ element, the part that explains why decisions are being made, as opposed to just a detailed implementation. This is equally important as reproducibility. And unfortunately, manuscript conventions for conciseness oftentimes push nitty gritty details out of the paper. With no regulated medium for these details to be archived, they more often than not end up cast into the wind and ultimately lost.

With this in mind, this document aims to explain both the _how_ and the _why_ of our research, by intermingling step-by-step reproducibility with in-depth commentary. This is because in our experience readers are much more interested in drawing inspiration for their own research than they in replicating our results, and this document is structured to reflect this.

### Goals

As stated, this document is a co-mingling of step-by-step reproducibility, and in-depth commentary. In my opinion (some of this is subjective), the document satisfies all of the following goals:

1. Step-by-step command-line instructions are provided for generating all raw, intermediate, and final data used in the paper, including all figures, tables, and numbers in the main and supplemental text.
2. There exists many 'checkpoint files', so that one may download intermediate data if they can't or don't want to follow the entirety of the workflow. This is especially important for steps requiring significant computational resources. For example, this workflow involves downloading roughly 4Tb of metagenomic and metatranscriptomic sequence data for a read recruitment experiment. Yet those unable or uninterested in allotting this amount of storage space can instead download the anvi'o databases (a few Gb) that result from this read recruitment experiment and continue with the more interesting aspects of the analysis.
3. I've avoided common pitfalls of replicability by paying extra attention to **(a)** creating isolated installation environments, **(b)** meticulous control over program versioning, **(c)** cross-platform compatibility, and **(d)** prohibited use of environmental variables that are specific to my personal environment.
4. Let's face it: there isn't enough space in the Methods section for all relevant details. Beyond what's required for reproducibility, this document houses information about methods that didn't make the cut, but are important nonetheless.
5. In my opinion, this document sets a standard for what transparency in computational science means, and why that's so valuable. This analysis is completely laid bare, so that it may be openly challenged, criticized, riffed on, and hopefully improved upon in the process. It likely contains errors--but that's okay. I didn't make this document to prove everything I did was right. I made this document to be educational.

### Summary

What can you expect to find in this document? Perhaps you're here for something specific. The following provides a summary of the analyses performed in this document, so you may jump around on a whim.

(FIXME insert a table of contents, like this one used for infant gut)

* [**Chapter I**: Genome-resolved Metagenomics](#chapter-i-genome-resolved-metagenomics)
* [**Chapter II**: Automatic Binning](#chapter-ii-automatic-binning)
* [**Chapter III**: Phylogenomics](#chapter-iii-phylogenomics)
* [**Chapter IV**: Pangenomics](#chapter-iv-pangenomics)
* [**Chapter V**: Metabolism Prediction](#chapter-v-metabolism-prediction)
* [**Chapter VI**: Microbial Population Genetics](#chapter-vi-microbial-population-genetics)
* [**Chapter VII**: Genes and genomes across metagenomes](#chapter-vii-genes-and-genomes-across-metagenomes)
* [**Chapter VIII**: From single-amino acid variants to protein structures](#chapter-viii-from-single-amino-acid-variants-to-protein-structures)

### Requirements

There are some skill and hardware requirements for following this workflow. There are also software requirements, however after completing the first step in this workflow, it is guaranteed that all of your software requirements are met.

#### Hardware

This study makes primary use of some large datasets. If you want to reproduce this study from scratch, you will need access to a high performance computing (HPC) cluster.

But if you don't have access to such a resource, don't worry. All of the computationally intensive steps in this study occur near the beginning of the workflow. Past Step FIXME, **all of the subsequent analysis is achievable on a laptop with FIXME Gb RAM and ~FIXME Gb storage space**.

I created several checkpoint data packs that you can download in order to jump into the remainder of the analysis without having completed earlier steps.

#### Skillset

This workflow will use language that assumes you have basic knowledge of UNIX (command-line) based scientific computing.

That means you know the basics of navigating to different files on your computer using commands like `cd` and `ls`. It means you shouldn't be surprised that if you have an R script called `test.R`, you could execute that script with the following syntax: `Rscript test.R`.

If you don't fall into this category of people, I would suggest checking out this great learning resource, [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/), which teaches you the basics of UNIX command-line computing.

If you are concerned whether your skills are up to snuff, don't worry. My suggestion would be to take a look at [Step FIXME: Setting up the required environment](FIXME). This is by far the hardest part of the workflow, so if you're able to complete this step and pass all the tests, the rest is child's play.

As you learned in the hardware requirements, recreating this study from scratch requires access to a high-performance computing cluster. If you plan to go this route instead of downloading some of the checkpoint data packs, you will need to know how to submit commands as jobs to your HPC cluster.

### Format

I tried to adopt a consistent "language" throughout the document, so please take the time to learn it.

#### Steps

FIXME I still need to figure out how I'm going to do this.

This workflow is composed of **steps** and **commands**. Each step is series of commands.

At the start of each step is a blue box that specifies the prerequisite steps and files needed to start the step. Alternatively, any checkpoint data packs will also be provided here, if one wants to jump into the section without having completed all previous sections. An example looks like this:

<div class="extra-info" markdown="1">
<span class="extra-info-header">Step X Info</span>
‣ **Prerequisite steps:** None FIXME add real example  
‣ **Checkpoint datapack:** None  
</div>

#### Commands

Each step contains commands that you should ideally be able to run in succession by copy-pasting them into your terminal prompt.

To distinguish these commands from all other commmentary and code, I have put each command inside a red (FIXME match description to final CSS) box. In addition to the command itself, each red box contains the time and storage requirements for running the command. If the command requires internet, I add the flag: '**Internet**: Yes'. If the command should be submitted as a job due to high memory/storage/time requirements, I add the flag: '**Cluster**: Yes'. Here is an example command:

<div class="extra-info" style="{{ command_style }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
anvi-split -C GENOMES \
           -c 03_CONTIGS/SAR11_clade-contigs.db \
           -p 06_MERGED/SAR11_clade/PROFILE.db \
           -o 07_SPLIT
```
‣ **Time:** 90 min  
‣ **Storage:** 28.3 Gb  
‣ **Cluster:** Yes (40 Gb)
</div>

If there are sections of commands that need to be edited manually by you, **required edits will be denoted within angular brackets (`<>`)**. For example,

<div class="extra-info" style="{{ command_style }}" markdown="1">
<span class="extra-info-header">Command #Y</span>

```bash
cd <A_DIRECTORY_YOU_LIKE>
mkdir kiefl_et_al_2021
```
‣ **Time:** <1 min  
‣ **Storage:** Minimal  
</div>

In this case, `<A_DIRECTORY_YOU_LIKE>` refers to a directory you decide. For example, I might replace it with `/Users/evan/Desktop`.

#### Help is a click away

**Purple keywords refer to anvi'o concepts that have help pages**. For example, {% include ARTIFACT name="contigs-db" %} refers to an anvi'o contigs database that can be clicked on to learn more about.

**Green keywords refer to anvi'o programs that have help pages**. For example, {% include PROGRAM name="anvi-gen-contigs-database" %} refers to an anvi'o program that can create an anvi'o contigs database.

The links within these documents lead to other anvi'o concepts, which lead to other anvi'o concepts. It's kind of like Wikipedia--you can get lost in it. Because many concepts already have rich documentation, I often forego in depth explanations, assuming you will click the link to learn more. 

## Directory

With that all out of the way, it is time for the workflow itself. Due to its size, I've broken it up into several parts, each one being a separate document. At the top of document you can find a directory for quickly switching between all of the different parts. Below is a description of each document.

### Chapter I: The prologue

This is what you're reading right now. In this section the philosophy, goals, expectations, requirements, and gameplan are all laid out.

[[Click here for Chapter I]]({{ site.url }}/data/anvio-structure/chapter-I)


### Chapter II: Configure your system

In this section, you'll set up your computational environment very precisely so that everything downstream runs as expected and you don't run into errors. This is the hardest chapter of the the workflow because exactness is key. But you'll be walked through each step, so don't worry.

[[Click here for Chapter II]]({{ site.url }}/data/anvio-structure/chapter-II)

### Chapter III: Build the data

In this section, you'll start with nothing but the internet, and build up all of the required data structures needed to reproduce all of the analyses in the paper.

[[Click here for Chapter III]]({{ site.url }}/data/anvio-structure/chapter-III)

### Chapter IV: Analyze the data

In this section, you'll utilize all of the data structures you've worked so hard to build in Chapter III, and reproduce each analysis in the study (or just the one's you care about).

[[Click here for Chapter IV]]({{ site.url }}/data/anvio-structure/chapter-IV)

### Chapter V: Reproduce every number

In this section, you'll be provided the required scripts to reproduce every mentioned number in the main text and supplemental text of the study.

[[Click here for Chapter V]]({{ site.url }}/data/anvio-structure/chapter-V)

