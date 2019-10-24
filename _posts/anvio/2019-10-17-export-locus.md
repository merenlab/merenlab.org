---
layout: post
authors: [matt, evan, alon]
title: "Extract loci from genomes and metagenomes with anvi'o!"
excerpt: "How to locate and extract target genetic loci from larger genetic context."
modified: 2019-10-08
tags: []
categories: [anvio]
comments: true
redirect_from:
  - /scg-taxonomy/
---

{% include _toc.html %}

{:.warning}
This tutorial is for `v6` and later versions of anvi'o. You can identify which version you have on your computer by typing `anvi-self-test --version` in your terminal.

{% include _project-anvio-version.html %}

## Introduction

Some genetic analyses call for the comparison of specific genetic loci between genomes. For example, one may be interested in investigating evidence for adaptive evolution of the lac operon between different _E. coli_ strains. The first step to this analysis would be to extract the various lac operon from a collection of _E. coli_ genomes. 

To address this example and other genomic loci analyses alike, we present the tool [anvi-export-locus](/software/anvio/vignette/#anvi-export-locus)!

Briefly, `anvi-export-locus` cuts out loci using two approaches: `default-mode` or what we call `flank-mode`. In the `default-mode`, the tool locates a designated anchor gene, then cuts upstream and downstream based on user-defined input. Notice that what is "upstream" and what is "downstream" is determined according to the direction of the anchor gene, i.e., if the anchor gene is in the reverse direction, then "upstream" would mean genes that have higher gene callers ids, and vice versa. On the other hand, `flank-mode` finds designated genes that define the left and right boundaries of the target locus, then cuts in between them. Genes to locate locus anchors or flanking genes are defined through their specific ids in anvi'o or through `search-terms` that query functional annotations or HMM hits stored in your contigs database.

To explore the functionality of `anvi-export-locus`, let's extract the lac operon from the larger genomic context of _E. coli_ genomes!

<div class="extra-info" markdown="1">

<span class="extra-info-header">Using anvi-export-locus with external gene calls</span>
If you are using [external gene calls](https://github.com/merenlab/anvio/issues/368) when generating your contigs databases, then you should make sure genes that are next to each other on a contig, have sequenctial gene caller ids. If you are using [Prodigal](https://github.com/hyattpd/Prodigal), then you have nothing to worry, since this is how Prodigal behaves and this is how we tested `anvi-export-locus`. To clarify this point, the way `anvi-export-locus` finds the genes for it to cut is by using the gene callers id. For example, if you are `flank-mode` and the the upstream and downstream genes that you provided hit gene caller ids 40 and 50, respectively, then `anvi-export-locus` will return the eleven genes with ids between 40 and 50. Intuitively, you would expect these to be next to each other on the chromosome, but there is nothing stopping you from using external gene calls that don't abide this intuition, and anvi'o is not going to check whether that is the case or not.
</div>


{% include _join-anvio-slack.html %}

## Downloading _E. coli_ genomes

First, let's download Genbank files for a few representative _E. coli_ strains:

{:.notice}
Alon has a great tutorial [here](http://merenlab.org/2019/03/14/ncbi-genome-download-magic/) describing how to automagically download genomes from NCBI. We'll be using this tool to download a few _E. coli_ genomes. If you have any questions regarding downloading genomes, please refer to [Alon's tutorial](http://merenlab.org/2019/03/14/ncbi-genome-download-magic/).

```{bash}
# Set working directory variable for later
WD=$(pwd)

# Download Genbank files
ncbi-genome-download bacteria \
                     --assembly-level chromosome,complete \
                     --genus Escherichia \
                     --metadata metadata.txt \
                     --refseq-category reference
```

Next, we'll make FASTAs, external gene calls, and functional annotations for all the Genbanks we just downloaded:
```{bash}
anvi-script-process-genbank-metadata -m metadata.txt \
                                     --output-dir ecoli \
                                     --output-fasta-txt ecoli.txt
```

## Generate contigs DBs
Now we need to get the fasta files into an anvi'o friendly format. To do this, we'll convert the files we made in the previous step into contigs DBs using the [anvi'o contigs workflow](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#contigs-workflow). 

First, make a json file for the anvio workflow called `contigs.json`
```{bash}
{
    "fasta_txt": "ecoli.txt"
}
```

Then run the contigs workflow! 

{:.notice}
This step may take a while depending on your computational resources. If you have any questions about running anvi'o workflows please refer to this tutorial [here](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#a-general-introduction-to-essentials). If you access to an HPC or cluster computer, check out additional parameters  [here](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#running-workflows-on-a-cluster).
```{bash}
anvi-run-workflow -w contigs \
                  -c contigs.json \
                  --additional-params \
                            --jobs 6 \
                            --resources nodes=6
```

## Extract lac operon

We should now have contigs DBs for our genomes. 

Now that we have six representative _E. coli_ examples, we will use a combination of the `default-mode` and `flank-mode` to cut out the genomic neighborhood around the lac operon and then trim the contig to just contain the target operon. 

First, we will use `default-mode`. This requires the user to provide a `--search-term` and `--num-genes` parameter. The `--search-term` will act as an anchor gene to locate the locus within the contigs provided. In this case, we will use the lacZ gene to locate the lac operon in our  _E. coli_ genomes. Once the anchor gene is located `--num-genes X,Y` will instruct `anvi-export-locus` to cut `X` gene(s) upstream and `Y` gene(s) downstream of the designated anchor gene.

Let's get cutting!

## Default mode

First, we'll use `default-mode` to extract the general genomic neighborhood around the lac operon in each genome
```{bash}
mkdir 03_LOCI

cd 03_LOCI

for GENOME in `ls "${WD}"/02_CONTIGS/*-contigs.db`;
do
    FNAME=$(basename "${GENOME}" -contigs.db)
    anvi-export-locus -c "${GENOME}" \
                      --num-genes 10,10 \
                      --search-term "lacZ" \
                      -O "${FNAME}"_lac_locus;
done
```


`--search-term` is NOT case sensitive unless you surround your term in quotes (e.g. `--search-term "lacZ"`)

Here is a visual representation of how `anvi-export-locus` found the anchor gene "lacZ" the cuts 10 genes upstread and downstream.
[![export-locus-defaultmode](/images/export-locus-defaultmode.png)](export-locus-defaultmode.png){:.center-img .width-100}

## Flank-mode

Awesome, now we have some smaller contigs that contain the lac operon. BUT, we also grabbed some extra genes that don't belong to the operon. Let's use `--flank-mode` to trim the loci to just contain the lac operon.

To do this, give `anvi-export-locus` two flanking `--search-term`s: lacI and lacA
```{bash}
for GENOME in `ls "${WD}"/03_LOCI/*.db`;
do
    FNAME=$(basename "${GENOME}" _lac_locus_0001.db)
    anvi-export-locus -c "${GENOME}" \
                      --flank-mode \
                      --search-term "lacI","lacA" \
                      -O "${FNAME}"_lac_locus_clean;
done
```

{:.notice}
`--flank-mode` requires flanking genes to be single copies in the contig it's searching. If your locus of interest does not have fixed coordinates in your genomes or metagenomes, you may need to adjust the `-search-term`s on a case by case bases. 

Here is a visual representation of how `flank-mode` cuts out a locus using flanking genes.
[![export-locus-defaultmode](/images/export-locus-flankmode.png)](export-locus-flankmode.png){:.center-img .width-100}

## Conclusion

[anvi-export-locus](/software/anvio/vignette/#anvi-export-locus) is a flexible tool that allows you to extract genomic loci from genomes and metagenomes. In this tutorial, we looked at the classic lac operon in _E. coli_ genomes, but this tool can also be unleashed on any genomic loci and even metagenomic assemblies! For instance, one could extract specific cellulose synthesis operons from soil metagenomes.

We hope you find amazing applications for this tool. If you have a suggestion or question please do not hesitate to contact us. Also, please report any bugs as an issue on the anvi'o [Github repository](https://github.com/merenlab/anvio).

{% include _join-anvio-slack.html %}

