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
This post also doubles as a reproducible workflow. Feel free to download the associated datapack at [FIXME INSERT DATAPACK LINK]() and follow along with the commands. (Or, go your own way and explore the data yourself!) The commands below were written for anvi'o `v7.1`. If you have a newer version of anvi'o and find a command that is not working - sorry. We don't always retroactively update these posts as anvi'o evolves. Please send us a message and we'll see if we can fix it.  

[TODO: change docs links to v7.1 instead of main]

## A focus on nitrogen fixation

About one year ago, I was looking for nitrogen fixation genes in ocean metagenomes. To tell you the truth, it wasn't so much that I cared a lot about nitrogen fixation back then (I had only just started getting into the whole metabolism thing), and more that I really just wanted to test out [this new program that I wrote](https://merenlab.org/software/anvio/help/main/programs/anvi-estimate-metabolism/) in as many different metagenome types as I could. So when my much-more-knowledgable colleagues suggested that I use this tool to search for nitrogen-fixing populations in a non-conventional environment - polar ocean metagenomes - I said, "Why not?"

For anyone who is as clueless as I was about nitrogen fixation, here is a very light summary of the background. Nitrogen fixation is the process of converting gaseous nitrogen (N<sub>2</sub>) into the more biologically-usable form ammonia (NH<sub>3</sub>). The resulting ammonia can then be further converted into other bioavailable compounds like nitrate and nitrite. Since nitrogen is an essential component of many biological molecules (amino acids, anyone?), nitrogen fixation is a fairly important process that supports life, in general. Only some ocean microbes, called marine diazotrophs, have the genes that allow them to fix nitrogen, which are (usually) encoded in the _nif_ operon. ([Sohm et al 2011](https://www.nature.com/articles/nrmicro2594))

Why are polar oceans a non-conventional place to find nitrogen fixers? Well, the majority of nitrogen fixing microbes have been found in non-polar oceans, perhaps because the diazotrophic cyanobacteria that are typically studied tend to be found in warmer waters ([Stal 2009](https://doi.org/10.1111/j.1758-2229.2009.00016.x)), or perhaps because the polar oceans are just not as well-studied as the other oceans. Regardless, in recent times, there have been several reports of nitrogen fixation happening in the Arctic and Antarctic Oceans [FIXME CITATIONS - https://www.frontiersin.org/articles/10.3389/fmicb.2020.596426/full, https://www.nature.com/articles/s41561-020-00651-7, https://www.pnas.org/content/115/52/13371]

To get back to the story, in early 2020 I started looking for nitrogen fixation genes in some polar ocean metagenomes, graciously shared by our collaborators at the [IUEM-Brest](https://www-iuem.univ-brest.fr/?lang=en) (who will hopefully share their findings and dataset in an exciting study in the future). And I found nothing. Nada. Not a single nitrogen fixation gene. I promise that I tried very hard, but you cannot find things if they are not there.

Thankfully, Meren soon brought to my attention a recently published dataset of Arctic and Antarctic ocean metagenomes by [Cao et al](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00826-9). These brave scientists faced the cold to bring the marine science community 60 new samples from 28 different locations in the polar oceans.  Their comparative analyses demonstrated that polar ocean microbial communities are distinct from non-polar ones, both in terms of their taxonomic diversity and their gene content. They also reconstructed 214 metagenome-assembled genomes, 32 of which were enriched in the polar oceans according to read recruitment analyses - the authors analyzed some metabolic pathways in these MAGs, but notably did not check for nitrogen fixation. So here was the opportunity for me and `anvi-estimate-metabolism` to try again.

Given my complete lack of success at finding anything nitrogen-fixy in the previous dataset, I was rather skeptical that this new dataset, while exciting, would turn out any different. But I gave it a shot, and - **spoilers - it worked! I not only found complete nitrogen fixation pathways in 4 individual Arctic Ocean metagenomes, but I also used that information to perform targeted binning of a microbial population genome containing one such pathway. It was a very happy moment, and I am excited to share with you how I did it.

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
    -O $name \
    --kegg-output-modes modules,kofam_hits; \
done < <(tail -n+2 metagenomes.txt)
```
{:.notice}
If you are determined to run this loop yourself, it is probably only possible on a high-performance computing cluster, in which case you will almost certainly have to modify the command to give each `anvi-estimate-metabolism` job more memory.

What this loop does is read each line of the `metagenomes.txt` file, except for the first one (the `tail -n+2` command skips the first line). Each non-header line in the file contains the name of the metagenome sample (which gets placed into the `$name` variable) and the path to its contigs database (which gets placed into the `$path` variable). Therefore, `anvi-estimate-metabolism` gets run on each contigs database in metagenome mode, and the resulting output files (two per sample) are prefixed with the sample name.

It _is_ possible to run `anvi-estimate-metabolism` on more than one contigs database at a time, using multi-mode, which you can read about on the %{anvi-estimate-metabolism}s help page. However, I did not do this here because I wanted the output for each sample to be printed to a separate output file, for purely organizational purposes.

You can download the resulting output files from [this link](FIXME FIGSHARE LINK). Once you do that, you will notice that there are 32 text files, two for each metagenome assembly. Let's take a look at the first few lines of the `modules` file for sample N02:
```bash
head -n 4 N02_modules.txt
```

You should see something like this:

unique_id | contig_name | kegg_module | module_name | module_class | module_category | module_subcategory | module_definition | module_completeness | module_is_complete | kofam_hits_in_module | gene_caller_ids_in_module | warnings
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
0 | c_000000008738 | M00546 | Purine degradation, xanthine => urea | Pathway modules | Nucleotide metabolism | Purine metabolism | "(K00106,K00087+K13479+K13480,K13481+K13482,K11177+K11178+K13483) (K00365,K16838,K16839,K22879) (K13484,K07127 (K13485,K16838,K16840)) (K01466,K16842) K01477" | 0.5 | False | K13485,K16842,K07127 | 121398,121397,121396 | None
1 | c_000000000052 | M00001 | Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "(K00844,K12407,K00845,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)" | 0.4 | False | K00134,K00873,K01624,K00927 | 8515,8519,8523,8454 | None
2 | c_000000000052 | M00002 | Glycolysis, core module involving three-carbon compounds | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)" | 0.5 | False | K00134,K00873,K00927 | 8515,8519,8454 | None

This is a [modules mode](https://merenlab.org/software/anvio/help/main/artifacts/kegg-metabolism/#modules-mode) output file from `anvi-estimate-metabolism` (which is the default output type). Since we ran the program in metagenome mode, each row of the file describes the completeness of a metabolic module within one contig of the metagenome. What this means is that every KOfam hit belonging to this pathway (listed in the `kofam_hits_in_module` column) was present on the _same contig_ in the metagenome assembly. This is important, because metagenomes contain the DNA sequences of multiple organisms, so the only time that we can be sure two genes go together within the same population genome is when they are assembled together onto the same contig sequence.

If right now you are thinking, "But wait... if we only focus on the genes within the same contig, many metabolic pathways will have completeness scores that are too low," then you are exactly correct. It is likely that most metabolic pathways from the same genome will be split across multiple contigs, so they will end up in different lines of this file. In the example above, contig `c_000000008738` contains 50% of the KOs required for the purine degradation module `M00546`, but perhaps the other KOs in the pathway (such as `K01477`) are also belonging to whatever microbial population this is, just on a different contig. Putting many contigs together to match up the different parts of the pathway, while making sure that you are not producing a chimeric population, is a task that would require careful binning.

Luckily for us, the nitrogen fixation module from KEGG has a couple of helpful characteristics. First, it contains only 3 genes, and second, it is encoded in an operon, so those genes are located close together in any given genome sequence. These two things make it much more likely that the entire module will end up within a single contig in our metagenome assemblies, which means it will be relatively easier to find a complete nitrogen fixation module in our metabolism estimation output files.

But before we dive into our search, let's quickly discuss what the nitrogen fixation operon, and its corresponding KEGG module, look like.

## Nitrogen fixation - KEGG vs reality

The KEGG module for nitrogen fixation is [M00175](https://www.genome.jp/module/M00175), and it looks like this:

[TODO: Insert M00175 screenshot]

You need a nitrogenase enzyme complex to convert nitrogen gas to ammonia, and there are only two possible versions of this complex: the "molybdenum-dependent nitrogenase" protein complex encoded by genes _nifH_ (K02588), _nifD_ (K02586), and _nifK_ (K02591) of the _nif_ operon; or the "vanadium-dependent nitrogenase" protein complex encoded by genes _vnfD_ (K22896), _vnfK_ (K22897), _vnfG_ (K22898), and _vnfH_ (K22899), of the (you guessed it) _vnf_ operon. The latter complex has been isolated from soil bacteria and is known to be an alternative nitrogenase that is expressed when molybdenum is not available ([Lee et al 2009](https://doi.org/10.1073/pnas.0904408106), [Bishop et al 1980](https://doi.org/10.1073/pnas.77.12.7342)). We're just going to ignore it, because I have yet to see it in any ocean samples.

So that means we effectively care about only _nifHDK_ in this module. But wait. While _nifHDK_ represent the catalytic components of the nitrogenase enzyme, it turns out that there are a few other genes required to produce the essential FeMo-cofactor and incorporate it into this complex. At minimum, the extra genes required are _nifE_ ([K02587](https://www.genome.jp/dbget-bin/www_bget?ko:K02587)), _nifN_([K02592](https://www.genome.jp/dbget-bin/www_bget?ko:K02592)), and _nifB_ ([K02585](https://www.genome.jp/dbget-bin/www_bget?ko:K02585))([Dos Santos et al 2012](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-162)). That means we need to find six genes - preferably on the same contig - within a metagenome sample in order to be confident that the metagenome includes a nitrogen-fixing population.

Just so you have a picture of what this should look like, here is a diagram of the _nif_ operon in the _Azotobacter vinelandii_ genome sequence, from [Dos Santos et al 2012](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-162)).

[TODO: insert Fig 1 from https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-162)]

The catalytic genes - all the ones from module M00175 - are located next to each other on the bacterial chromosome. The other required biosynthetic genes are located farther along, with _nifE_ and _nifN_ expressed under the same promoter and _nifB_ isolated from the rest of the genes and expressed under its own promoter.

We expect to see this general pattern reflected in the Arctic Ocean metagenome assemblies, meaning that genes _nifHDKEN_ are the most likely to end up on the same contig. If you keep reading, you will see that this is indeed the case!

Without further ado, let's take a look at the data.

## Looking for evidence of a nitrogen-fixing population

Though M00175 only contains the catalytic portion of our required _nif_ gene set, it is a good starting point for our search. If we look for this module in our metabolism estimation results, we can find out which contig(s) it is located on and use that to guide our search for the remaining genes.

You can use the following `bash` code to search for lines describing M00175 in all metabolism estimation `modules` mode outputs. The code filters the output so that it contains only those lines which have a score of 1.0 in the `module_completeness` column, meaning that all 3 _nifHDK_ genes are located on the same contig. It further filters the output to contain only the columns describing 1) the file name and line of file where M00175 was found, 2) the contig name, 9) the completeness score, 11) the list of KO hits that we found from this module, and 12) the corresponding gene caller IDs of these hits.

```bash
grep M00175 *_modules.txt | awk -F'\t' '$9 == 1.0' | cut -f 1,2,9,11,12
```

Your output should look like this:

N06-contigs_modules.txt:45627 | c_000000000414 | 1.0 | K02586,K02591,K02588 | 35075,35076,35074
N07-contigs_modules.txt:7413 | c_000000004049 | 1.0 | K02586,K02591,K02588 | 94224,94225,94223
N07-contigs_modules.txt:31467 | c_000000000073 | 1.0 | K02586,K02591,K02588 | 14638,14637,14639
N22-contigs_modules.txt:44057 | c_000000000122 | 1.0 | K02586,K02591,K02588 | 16856,16857,16855
N25-contigs_modules.txt:11798 | c_000000000104 | 1.0 | K02586,K02591,K02588 | 13919,13920,13918

These are promising results! The complete M00175 module was found in 4 different Arctic Ocean samples (there are two different instances in sample N07).

I encourage you to look through the other instances of this module in the output files. If you do this, you will see that some metagenomes appear to have all three of these genes split across multiple contigs (could they be contigs from the same genome?). For instance, here is a pair of contigs from sample N22:

N22-contigs_modules.txt:35879 | c_000000000861 | 0.3333333333333333 | K02588 | 43430
N22-contigs_modules.txt:49457 | c_000000003717 | 0.6666666666666666 | K02591,K02586 | 84130,84129

_nifH_ is on one contig and _nifDK_ are on the other. I think it is likely that these two contigs go together, because it seems unlikely that a genome would have one of these genes from this operon and not the rest (though it could happen, of course. Things like prophages and transposons often destroy our expectations for microbial genomes).

All in all, as you examine the estimation results for these 16 metagenomes, you should find that 9 of them have at least a partial copy of M00175, and 5 of those contain a complete set of _nifHDK_.

Of course, as we discussed earlier, there are 3 other genes that we need to find alongside _nifHDK_ in order to be sure that we have a microbial population capable of fixing nitrogen. KEGG may not have put these genes in M00175, but it does have a KOfam profile for each one of _nifENB_ - those KOs are K02587, K02592, and K02585. To search for these, we turn to our `kofam_hits` mode output files.

We will focus on the five samples that contain _nifHDK_, which are N06, N07, N22, N25, and N38.

[TODO: write code and output examples for section on kofam hits]

In sample N06, we found a complete M00175 module on contig `c_000000000414`. From the `kofam_hits` output, we can see that _nifE_ and _nifN_ are also on that contig, while _nifB_ is on a separate one (contig `c_000000001170`). This arrangement makes sense, since _nifB_ is the gene expected to be farthest away from the others based on the _A. vinelandii_ genome we looked at earlier. And since all six of the required _nif_ genes are present, it seems very likely that this metagenome contains a legitimate nitrogen-fixing population!

```
FIXME
```

In sample N07, we find all six genes on the same contig, `c_000000000073`. This is even better. There is a nitrogen-fixing population here for sure (and there may even be two different ones, considering that contig `c_000000004049` also contains a complete M00175 and there is a second set of the _nifENB_ genes spread across three different contigs).

```
FIXME
```

Sample N22 has a similar situation to N06, with 5 out of 6 genes appearing on contig `c_000000000122`. It also appears to have a second set of these genes spread across multiple contigs, just as in N07.

```
FIXME
```

You can take a look at N25 and N38 yourself. N25 should have at least one copy of all six genes (and 5/6 on the same contig), but N38 should be missing _nifN_.

At this point, we can be fairly confident that there are nitrogen-fixing populations in samples N06, N07, N22, and N25. The natural question to ask next is - what are they?

## Determining population identities

Since we are working with individual contigs and not full genomes right now, a good strategy to figure out what these populations could be is to use [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) to see if there is anything similar to these contigs in the NCBI database.

Well, because I am nice, I extracted the relevant contig sequences from these metagenome assemblies for you. You will find them in the [TODO MAKE FILE] file in the datapack.