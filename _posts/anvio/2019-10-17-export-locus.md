---
layout: post
authors: [matt, meren]
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

{: .notice}
This tutorial is for `v6` and later versions of anvi'o. You can identify which version you have on your computer by typing `anvi-self-test --version` in your terminal.

## Introduction
Some genetic analyses call for the comparison of specific genetic loci between genomes. For example, one may be interested in investigating evidence for adaptive evolution of the lac operon between different E. coli strains. The first step to this analysis would be to extract the various lac operon from a collection of E. Coli genomes. 

Today, we'll use the tool `anvi-export-locus` to extract the lac operon from the larger genomic context of E. coli genomes!

Briefly, `anvi-export-locus` cuts out loci using two approaches: `default-mode` or what we call `flank-mode`. In the `default-mode`, the tool locates a designated anchor gene, then cuts upstream and downstream based on user-defined input. `Flank-mode`, on the other hand, finds designated genes that surround the target locus, then cuts in between them. Genes to locate locus anchors or flanking genes are defined through their specific ids in anvi'o or through `search-terms` that query functional annotations or HMM hits stored in your contigs database!

Let's get started.

## Download E. Coli genomes

Alon has a great tutorial [here](http://merenlab.org/2019/03/14/ncbi-genome-download-magic/) describing how to automagically download genomes from NCBI. We'll be using this to download a few E. Coli genomes today. If you have any questions regarding downloading genomes, please refer to [Alon's tutorial](http://merenlab.org/2019/03/14/ncbi-genome-download-magic/).

First, we'll download Genbank files for a few representative E. Coli strains.
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

Next, we'll make FASTAs, external gene calls, and functional annotations for all the Genbanks we just downloaded.
```{bash}
anvi-script-process-genbank-metadata -m metadata.txt \
                                     --output-dir ecoli \
                                     --output-fasta-txt ecoli.txt
```

## Generate contigs DBs
Now let's run the contigs workflow to create contigs DBs for each of our genomes.

First, we'll make a json file for the anvio workflow called `contigs.json`
```{bash}
{
    "fasta_txt": "ecoli.txt"
}
```

Then run the workflow! This step may take a while depending on your computational resources. If you have any questions about running anvi'o workflows please refer to [here](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/).
```{bash}
anvi-run-workflow -w contigs \
                  -c contigs.json \
                  --additional-params \
                  --cluster \
                        "clusterize \
                            -j={rule} \
                            -o={log} \
                            -e={log} \
                            -n={threads} \
                            -x" \
                   --jobs 10 \
                   --resource nodes=40
```

## Extract lac operon
We should now have contigs DBs for our E. Coli genomes. 

Today with our six representative E. Coli examples we will use a combination of the `default-mode` and `flank-mode` to cut out the genomic neighborhood around the lac operon and then trim the contig to just contain the target operon. 

First, we will use `default-mode`. This requires you to provide a `-search-term` and `--num-genes` parameter. The `-search-term` will act as an anchor gene to locate the locus within the contigs provided. In this case, we will use the lacZ gene (β-galactosidase which cleaves lactose) to locate the lac operon in our  E. coli genomes. Once the anchor gene is located `--num-genes X,Y` will instruct `anvi-export-locus` to cut `X` gene(s) upstream and `Y` gene(s) downstream of the designated anchor gene.

Let's get cutting!

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

**NOTE**
`--search-term` is NOT case sensitive unless you surround your term in quotes (e.g. `--search-term "lacZ"`)

Awesome we have some smaller contigs that contain the lac operon, BUT we also grabbed some extra genes that don't belong to the operon. Let's use `--flank-mode` to trim the loci to just contain the lac operon.

To do this, we will give `anvi-export-locus` two flanking `--search-terms': mhpR and lacA
```{bash}
for GENOME in `ls "${WD}"/03_LOCI/*.db`;
do
    FNAME=$(basename "${GENOME}" _lac_locus_0001.db)
    anvi-export-locus -c "${GENOME}" \
                      --flank-mode \
                      --search-term mhpR,"lacA" \
                      -O "${FNAME}"_lac_locus_clean;
done
```

**NOTE**
`--flank-mode` requires flanking genes to be single copies in the contig it's searching. If your locus of interest does not have fixed coordinates in your genomes or metagenomes, you may need to adjust the `-search-terms` on a case by case bases. 

## Conclusion

`anvi-export-locus` is a flexible tool that allows you to extract genomic loci from genomes and metagenomes. Today we looked at the classic lac operon in E. Coli genomes, but this tool can also be unleashed on metagenomic assemblies! For instance, one could extract specific cellulose synthesis operons from soil metagenomes.

We hope you find amazing applications for this tool. If you have a suggestion or question please do not hesitate to contact us! Also, please report any bugs as an issue on the anvi'o [Github repository](https://github.com/merenlab/anvio).
