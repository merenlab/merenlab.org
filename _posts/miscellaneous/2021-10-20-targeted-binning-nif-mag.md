---
layout: post
authors: [iva]
title: "Targeted binning of a novel nitrogen-fixing population from the Arctic Ocean"
excerpt: "Using metabolism estimation to go fishing in some ocean metagenomes."
modified: 2021-10-20
tags: [metabolism]
categories: [miscellaneous, anvio]
comments: true
redirect_from:
  - /targeted-binning-nif-mag/
---

Public sequencing datasets are beautiful things. They open a world of endless possibilities to scientists everywhere simply by being downloadable. It doesn't matter what they were originally sequenced for, because curious people can always find something to explore in them. Stool metagenomes from microbiome studies like [this one](https://doi.org/10.1016/j.cub.2015.04.055) might eventually be used [to find human gut-associated plasmids](https://www.biorxiv.org/content/10.1101/2020.11.01.361691v1.full). Genomes and MAGs from sequencing projects all over the globe come together in the phylogenies and pangenomes of countless papers. Today's [time-series](https://pubmed.ncbi.nlm.nih.gov/22936250/) can turn into tomorrow's [tutorial](https://merenlab.org/tutorials/infant-gut/). That's the magic of open science.

At the same time, the curse of open science is that sequencing projects designed and carried out by other people are not often obviously or directly usable in your own projects. Usually, you don't find what you are looking for even if you have the right tools to look for it. But sometimes, you get lucky. This post is about one of those times.

If you are interested in learning about how to leverage anvi'o's metabolism estimation capabilities to go fishing through your data, or if you just get a kick out of success stories in the realm of open science and data reusability, keep on reading.

{:.notice}
This post also doubles as a reproducible workflow. Feel free to download the associated datapack at [FIXME INSERT DATAPACK LINK]() and follow along with the commands. (Or, go your own way and explore the data yourself!) The commands below were written for anvi'o `v7`. If you have a newer version of anvi'o and find a command that is not working - sorry. We don't always retroactively update these posts as anvi'o evolves. Please send us a message and we'll see if we can fix it.  

(TODO: check v7 compatibility of commands)

## A focus on nitrogen fixation

About one year ago, I was looking for nitrogen fixation genes in ocean metagenomes. To tell you the truth, it wasn't so much that I cared a lot about nitrogen fixation back then (I had only just started getting into the whole metabolism thing), and more that I really just wanted to test out [this new program that I wrote](https://merenlab.org/software/anvio/help/7/programs/anvi-estimate-metabolism/) in as many different metagenome types as I could. So when my much-more-knowledgable colleagues suggested that I use this tool to search for nitrogen-fixing populations in a non-conventional environment - polar ocean metagenomes - I said, "Why not?"

For anyone who is as clueless as I was about nitrogen fixation, here is a very light summary of the background. Nitrogen fixation is the process of converting gaseous nitrogen (N<sub>2</sub>) into the more biologically-usable form ammonia (NH<sub>3</sub>). The resulting ammonia can then be further converted into other bioavailable compounds like nitrate and nitrite. Since nitrogen is an essential component of many biological molecules (amino acids, anyone?), nitrogen fixation is a fairly important process that supports life, in general. Only some ocean microbes, called marine diazotrophs, have the genes that allow them to fix nitrogen, which are encoded in the _nif_ operon. ([Sohm et al 2011](https://www.nature.com/articles/nrmicro2594))

Why are polar oceans a non-conventional place to find nitrogen fixers? Well, the majority of nitrogen fixing microbes have been found in non-polar oceans, perhaps because the diazotrophic cyanobacteria that are typically studied tend to be found in warmer waters ([Stal 2009](https://doi.org/10.1111/j.1758-2229.2009.00016.x)), or perhaps because the polar oceans are just not as well-studied as the other oceans. Regardless, in recent times, there have been several reports of nitrogen fixation happening in the Arctic and Antarctic Oceans [FIXME CITATIONS - https://www.frontiersin.org/articles/10.3389/fmicb.2020.596426/full, https://www.nature.com/articles/s41561-020-00651-7, https://www.pnas.org/content/115/52/13371]

To get back to the story, in early 2020 I started looking for nitrogen fixation genes in some private polar ocean metagenomes, graciously shared by our collaborators at FIXME (who will hopefully share their findings in an exciting study in the future). And I found nothing. Nada. Not a single nitrogen fixation gene. I promise that I tried very hard, but you cannot find things if they are not there.

Thankfully, Meren soon brought to my attention a recently published dataset of Arctic and Antarctic ocean metagenomes by [Cao et al](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00826-9). These brave scientists faced the cold to bring the marine science community 60 new samples from 28 different locations in the polar oceans.  Their comparative analyses demonstrated that polar ocean microbial communities are distinct from non-polar ones, both in terms of their taxonomic diversity and their gene content. They also reconstructed 214 metagenome-assembled genomes, 32 of which were enriched in the polar oceans according to read recruitment analyses - the authors analyzed some metabolic pathways in these MAGs, but notably did not check for nitrogen fixation. So here was the opportunity for me and `anvi-estimate-metabolism` to try again.

Given my complete lack of success at finding anything nitrogen-fixy in the previous dataset, I was rather skeptical that this new dataset, while exciting and new, would turn out any different. But I gave it a shot, and - **spoilers - it worked! I not only found complete nitrogen fixation pathways in 4 individual Arctic Ocean metagenomes, but I also used that information to perform targeted binning of a microbial population genome containing one such pathway. It was a very happy moment, and I am excited to share with you how I did it.

Let's go through this analysis together :)

## Estimating metabolism in Arctic Ocean metagenomes

To start, we need metagenome assemblies of the Arctic Ocean samples from Cao et al's dataset. I am fortunate to be colleagues with [Matt Schechter](https://orcid.org/0000-0002-8435-3203), an awesome microbiologist who knows way more about oceans than I do, and who also happened to be interested in this dataset. He downloaded the samples and made single assemblies of them using the software [IDBA-UD](FIXME CITATION) as part of [the anvi'o metagenomic workflow](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#metagenomics-workflow). We are all benefiting from his hard work today - thanks, Matt!

We won't look at all 60 samples from the Cao et al paper, only 16 of their surface Arctic Ocean samples [FIXME CHECK IF SURFACE].

{:.notice}
You do not need the metagenome assemblies to follow the rest of this blog post, but you would like to have access to the 16 assemblies I am talking about, you can download their contigs databases [here](FIXME FIGSHARE LINK). But be warned - they will take up 6.5 GB of space on your computer.

{:.warning}
When we downloaded these samples, we assigned different (shorter) names to them, so the sample names I will discuss below are different from the ones in the Cao et al paper. If you want to know the correspondence between our sample names and those in the paper, check out the `FIXME` file in the datapack. [TODO: create sample correspondence file]

The first thing that I did with those 16 assemblies was run %{anvi-estimate-metabolism}s in metagenome mode. I will show you the commands that I used to do this, but I won't ask you to do it yourself, because it takes quite a long time (and currently requires an obscene amount of memory, for which I deeply apologize). I created a %{metagenomes}s file containing the names and %{contigs-db}s paths of each sample called `metagenomes.txt`, and I wrote a bash loop to estimate metabolism individually on each sample:

```bash
while read name path; \
do \
    anvi-estimate-metabolism -c $path \
    --metagenome-mode \
    -O $name; \
done < <(tail -n+2 metagenomes.txt)
```
{:.notice}
If you are determined to run this loop yourself, it is probably only possible on a high-performance computing cluster, in which case you will almost certainly have to modify the command to give each `anvi-estimate-metabolism` job more memory.

What this loop does is read each line of the `metagenomes.txt` file, except for the first one (the `tail -n+2` command skips the first line). Each non-header line in the file contains the name of the metagenome sample (which gets placed into the `$name` variable) and the path to its contigs database (which gets placed into the `$path` variable). Therefore, `anvi-estimate-metabolism` gets run on each contigs database in metagenome mode, and the resulting output file is prefixed with the sample name.

It _is_ possible to run `anvi-estimate-metabolism` on more than one contigs database at a time, using multi-mode, which you can read about on the %{anvi-estimate-metabolism}s help page. However, I did not do this here because I wanted the output for each sample to be printed to a separate output file, for purely organizational purposes.

You can download the resulting output files from [this link](FIXME FIGSHARE LINK). Once you do that, you will notice that there are 16 text files, one for each metagenome assembly. Let's take a look at the first few lines of the file for sample N02:
```bash
head -n 4 N02-contigs_modules.txt
```

You should see something like this:

unique_id | contig_name | kegg_module | module_name | module_class | module_category | module_subcategory | module_definition | module_completeness | module_is_complete | kofam_hits_in_module | gene_caller_ids_in_module | warnings
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
0 | c_000000008738 | M00546 | Purine degradation, xanthine => urea | Pathway modules | Nucleotide metabolism | Purine metabolism | "(K00106,K00087+K13479+K13480,K13481+K13482,K11177+K11178+K13483) (K00365,K16838,K16839,K22879) (K13484,K07127 (K13485,K16838,K16840)) (K01466,K16842) K01477" | 0.5 | False | K13485,K16842,K07127 | 121398,121397,121396 | None
1 | c_000000000052 | M00001 | Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "(K00844,K12407,K00845,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)" | 0.4 | False | K00134,K00873,K01624,K00927 | 8515,8519,8523,8454 | None
2 | c_000000000052 | M00002 | Glycolysis, core module involving three-carbon compounds | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)" | 0.5 | False | K00134,K00873,K00927 | 8515,8519,8454 | None

This is a [modules mode](https://merenlab.org/software/anvio/help/main/artifacts/kegg-metabolism/#modules-mode) output file from `anvi-estimate-metabolism` (which is the default output type). Since we ran the program in metagenome mode, each row of the file describes the completeness of a metabolic module within one contig of the metagenome. What this means is that every KOfam hit belonging to this pathway (listed in the `kofam_hits_in_module` column) was present on the _same contig_ in the metagenome assembly. This is important, because metagenomes contain the DNA sequences of multiple organisms, so the only time that we can be sure two genes go together within the same population genome is when they are assembled together onto the same contig sequence.

If right now you are thinking, "But wait... if we only focus on the genes within the same contig, many metabolic pathways will have completeness scores that are too low," then you are exactly correct. It is likely that most metabolic pathways from the same population genome will be split across multiple contigs, so they will end up in different lines of this file. In the example above, contig `c_000000008738` contains 50% of the KOs required for the purine degradation module `M00546`, but perhaps the other KOs in the pathway (such as `K01477`) are also belonging to whatever population this is, just on a different contig. Putting many contigs together to match up the different parts of the pathway, while making sure that you are not producing a chimeric population, is a task that would require careful binning.

