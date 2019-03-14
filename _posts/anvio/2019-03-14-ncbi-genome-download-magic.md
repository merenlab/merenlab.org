---
layout: post
title: "Easily download and process genomes from the NCBI databases"
excerpt: "How to combine genomes from NCBI in your pangenomic/phylogenomic analysis"
modified: 2019-03-14
categories: [anvio]
comments: true
authors: [alon]
---

{% include _toc.html %}

{:.notice}
This tutorial was written using anvi'o `v5.4`.

## Introduction

Often we find ourselves seeking to compare our genome/s of interest to similar genomes that are available online.
Here I describe how I have been going about it recently.

In this short post I will cover:
1. How I download genomes from NCBI (using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)).
2. How I process these genomes (using [anvi-script-process-genbank-metadata](http://merenlab.org/software/anvio/vignette/#anvi-script-process-genbank-metadata)).
3. Run the [contigs workflow](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#contigs-workflow) to generate a contigs database for each of these genomes (using `anvi-run-workflow`).

There has been a [recent post](http://merenlab.org/2018/12/01/combining-annotation-sources-for-pan/) covering this exact issue by Mike Lee. Here I just present a slightly different approach.

## Downloading the genomes

{:.notice}
Recently, when we want to compare one of our MAGs to similar things on NCBI, we have been loving [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download).
We also use it in combination with the helper script [gimme_taxa.py](https://github.com/kblin/ncbi-genome-download#contributed-scripts-gimme_taxapy).
Here I will show you how.

So imagine you generated a MAG, and taxonomic annotation places it in some phylum/class/order/whatever and you want to compare it to other genomes from this phylum/class/order/whatever from the NCBI, how do you go about it?

For example, let's assume we have a MAG that we resolve to the Candidate phylum Gracilibacteria (formerly [GN02](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1472358/)).
In order to use `ncbi-genome-download`, we need specific TaxIDs.
In order to find these we use the helper script.
First, we download the helper script:

```bash
wget https://raw.githubusercontent.com/kblin/ncbi-genome-download/master/contrib/gimme_taxa.py ./
```

We can now use this script to get the TaxID (along with other information):

```bash
python gimme_taxa.py Gracilibacteria \
                     -o GN02-TaxIDs-for-ngd.txt
```

Another equal way to do this, is to first go to NCBI's [Taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) and find out the taxon id for Gracilibacteria (by simply searching for the term "Gracilibacteria").
We find out that it is 363464.

Now we run the helper script to find the specific taxon id for each genome that belongs to taxa ID 363464:

```bash
python gimme_taxa.py 363464 \
                     -o GN02-TaxIDs-for-ngd.txt
```

This would give us an identical result.
The only reason I mention this, is that it sometimes might be better to use the taxa ID from NCBI to guarantee that you get what you want.

Let's look into the output:

```bash
$ head GN02-TaxIDs-for-ngd.txt
parent_taxid	descendent_taxid	descendent_name
363464	363504	uncultured Candidatus Gracilibacteria bacterium
363464	1130342	Gracilibacteria bacterium JGI 0000069-K10
363464	1130343	Gracilibacteria bacterium JGI 0000069-P22
363464	1151650	Gracilibacteria bacterium canine oral taxon 291
363464	1151651	Gracilibacteria bacterium canine oral taxon 323
363464	1151652	Gracilibacteria bacterium canine oral taxon 364
363464	1151653	Gracilibacteria bacterium canine oral taxon 394
363464	1226338	Gracilibacteria bacterium oral taxon 871
363464	1226339	Gracilibacteria bacterium oral taxon 872
```

In order to run ncbi-genome-download we just need a list of the speicific taxids (the above "descendent_taxid"), so we can run the helper script again like this:

```bash
python gimme_taxa.py Gracilibacteria \
                     -o GN02-TaxIDs-for-ngd-just-IDs.txt \
                     --just-taxids
```

We can now download all of these genomes like this:

```bash
ncbi-genome-download -t GN02-TaxIDs-for-ngd-just-IDs.txt \
                     bacteria \
                     -o GN02 \
                     -m GN02-NCBI-METADATA.txt \
                     -s genbank
```

{:.notice}
We specified `-s genbank` since there are no Gracilibacteria genomes on RefSeq, only on GenBank.
To learn more about the parameters, refer to the help menu of ncbi-genome-download:
`ncbi-genome-download -h`

{:.warning}
Pro tip: if you want to first see a list of what is going to be downloaded use the `-n` flag for a dry run.

## Process NCBI genomes and get thing ready for anvi-run-workflow

And now for the most beautiful part, once you are done downloading the files, a metadata file is also created by ncbi-genome-download.
Here is a glimpse to this metadata file:

```bash
$ column -t GN02-NCBI-METADATA.txt | head -n 3
assembly_accession  bioproject   biosample     wgs_master      excluded_from_refseq  refseq_category  relation_to_type_material  taxid            species_taxid  organism_name  infraspecific_name  isolate          version_status   assembly_level                            release_type                              genome_rep       seq_rel_date  asm_name     submitter    gbrs_paired_asm  paired_asm_comp  ftp_path     local_filename
GCA_000404985.1     PRJNA192474  SAMN02441616  ASNE00000000.1  derived               from             single                     cell             na             1130342        1130342             Gracilibacteria  bacterium        JGI                                       0000069-K10                               strain=JGI       0000069-K10   latest       Contig       Major            Full             2013/06/03   ASM40498v1       DOE          Joint                                                                               Genome                                                                               Institute                                                                            GCF_000404985.1                                                                      identical                                                                            ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/985/GCA_000404985.1_ASM40498v1    ./GN02/genbank/bacteria/GCA_000404985.1/GCA_000404985.1_ASM40498v1_genomic.gbff.gz
GCA_000404725.1     PRJNA192390  SAMN02440944  ASND00000000.1  derived               from             single                     cell             ;              partial        na                  1130343          1130343          Gracilibacteria                           bacterium                                 JGI              0000069-P22   strain=JGI   0000069-P22  latest           Contig           Major        Partial          2013/06/03   ASM40472v1                                                                          DOE                                                                                  Joint                                                                                Genome                                                                               Institute                                                                            na                                                                                   na                                                                                   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/725/GCA_000404725.1_ASM40472v1   ./GN02/genbank/bacteria/GCA_000404725.1/GCA_000404725.1_ASM40472v1_genomic.gbff.gz
```

We created a sctipt that takes this metadata file and generates a [fasta.txt](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt)
file in the format compatible with the anvi'o workflows. Here is how you run this script for our example:

```bash
anvi-script-process-genbank-metadata -m GN02-NCBI-METADATA.txt \
                                     -o GN02-NCBI-GENOMES \
                                     --output-fasta-txt GN02-fasta.txt
```

The output directory `GN02-NCBI-GENOMES` contains three files for each of the genomes that were downloaded:
1. FASTA files
2. External gene calls file - formated so that the [PGAP](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/) gene calls could be imported into an anvio contigs database.
3. External gene functions file - formatted properly to be imported into an anvi'o contigs database.

Here is what this folder looks like in our case:

```bash
$ ls GN02-NCBI-GENOMES/
Candidatus_Altimarinus_pacificus_JGI_0000068_E11_GCA_000405005.1-contigs.fa						Candidatus_Gracilibacteria_bacterium_GCA_002748005.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-contigs.fa			Candidatus_Gracilibacteria_bacterium_GCA_002748005.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-external-functions.txt		Candidatus_Gracilibacteria_bacterium_GCA_002749055.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-external-gene-calls.txt	Candidatus_Gracilibacteria_bacterium_GCA_002749055.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-contigs.fa		Candidatus_Gracilibacteria_bacterium_GCA_002749055.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-external-functions.txt	Candidatus_Gracilibacteria_bacterium_GCA_003242835.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-external-gene-calls.txt	Candidatus_Gracilibacteria_bacterium_GCA_003242835.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-contigs.fa			Candidatus_Gracilibacteria_bacterium_GCA_003242835.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-external-functions.txt		Candidatus_Gracilibacteria_bacterium_GCA_003251355.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-external-gene-calls.txt		Candidatus_Gracilibacteria_bacterium_GCA_003488655.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-contigs.fa						Candidatus_Gracilibacteria_bacterium_GCA_003488655.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-external-functions.txt				Candidatus_Gracilibacteria_bacterium_GCA_003488655.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-external-gene-calls.txt				Candidatus_Gracilibacteria_bacterium_GCA_003638805.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-contigs.fa						Candidatus_Gracilibacteria_bacterium_GCA_003638805.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-external-functions.txt				Candidatus_Gracilibacteria_bacterium_GCA_003638805.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-external-gene-calls.txt				Candidatus_Gracilibacteria_bacterium_GCA_003638815.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-contigs.fa				Candidatus_Gracilibacteria_bacterium_GCA_003638815.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-external-functions.txt		Candidatus_Gracilibacteria_bacterium_GCA_003638815.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-external-gene-calls.txt		Candidatus_Gracilibacteria_bacterium_GN02_872_GCA_003260325.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-contigs.fa				Candidatus_Gracilibacteria_bacterium_GN02_872_GCA_003260325.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-external-functions.txt		Candidatus_Gracilibacteria_bacterium_GN02_872_GCA_003260325.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-external-gene-calls.txt		Candidatus_Gracilibacteria_bacterium_GN02_873_GCA_003260345.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-contigs.fa								Candidatus_Gracilibacteria_bacterium_GN02_873_GCA_003260345.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-external-functions.txt						Candidatus_Gracilibacteria_bacterium_GN02_873_GCA_003260345.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-external-gene-calls.txt						Candidatus_Gracilibacteria_bacterium_HOT_871_GCA_002761215.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_GCA_002747975.1-contigs.fa								Candidatus_Gracilibacteria_bacterium_HOT_871_GCA_002761215.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_GCA_002747975.1-external-functions.txt						Candidatus_Gracilibacteria_bacterium_HOT_871_GCA_002761215.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_GCA_002747975.1-external-gene-calls.txt						Gracilibacteria_bacterium_JGI_0000069_K10_GCA_000404985.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_GCA_002748005.1-contigs.fa								Gracilibacteria_bacterium_JGI_0000069_P22_GCA_000404725.1-contigs.fa
```

The other output of `anvi-script-process-genbank-metadata` is the fasta.txt file. Here is how this file looks like in our case:

```bash
$ column -t GN02-fasta.txt | head
name                                                                                           path                                                                                                                                                     gene_functional_annotation                                                                                                                                           external_gene_calls
Candidatus_Altimarinus_pacificus_JGI_0000068_E11_GCA_000405005.1                               /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Altimarinus_pacificus_JGI_0000068_E11_GCA_000405005.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1       /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-contigs.fa       /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-external-functions.txt       /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1  /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-contigs.fa  /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-external-functions.txt  /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1        /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-contigs.fa        /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-external-functions.txt        /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1                             /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-contigs.fa                             /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-external-functions.txt                             /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1                              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-contigs.fa                              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-external-functions.txt                              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1            /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-contigs.fa            /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-external-functions.txt            /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-contigs.fa              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-external-functions.txt              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-external-gene-calls.txt
Candidatus_Gracilibacteria_bacterium_GCA_002746815.1                                           /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-contigs.fa                                           /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-external-functions.txt                                           /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-external-gene-calls.txt
```

We can now use the TXT file `GN02-fasta.txt` in a config file for a [contigs](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#contigs-workflow), [phylogenomics](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#phylogenomics-workflow) or [pangenomics](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#pangenomics-workflow) workflows (if you are not familiar with these workflows refer to the tutorials in the provided links <<).

### Excluding the gene calls from the fasta.txt file
When you specify an external gene calls file in your fasta.txt file, then `anvi-run-workflow` will use these gene calls instead of running prodigal.
In this case the "soucre" for the gene calls will be noted as "PGAP".
This means that we might end up having some of the genomes with gene calls from PGAP and others with gene calls from prodigal.
But because of the problem descirbed in this [github issue](https://github.com/merenlab/anvio/issues/565), when we run a pangenomic or phylogenomic analysis with anvi'o,
all the genomes must contain gene calls from the same source.

In order to bypass this issue, this is how I run `anvi-script-process-genbank-metadata`:

```bash
anvi-script-process-genbank-metadata -m GN02-NCBI-METADATA.txt \
                                     -o GN02-NCBI-GENOMES \
                                     --output-fasta-txt GN02-fasta.txt \
                                     --exclude-gene-calls-from-fasta-txt
```

The output directory would stay the same, but the TXT file `GN02-fasta.txt` now looks like this:

```bash
$ column -t GN02-fasta.txt | head
name                                                                                           path                                                                                                                                                     gene_functional_annotation
Candidatus_Altimarinus_pacificus_JGI_0000068_E11_GCA_000405005.1                               /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Altimarinus_pacificus_JGI_0000068_E11_GCA_000405005.1-contigs.fa
Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1       /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-contigs.fa       /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG12_big_fil_rev_8_21_14_0_65_38_15_GCA_002771035.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1  /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-contigs.fa  /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG17_big_fil_post_rev_8_21_14_2_50_48_13_GCA_002783265.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1        /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-contigs.fa        /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG18_big_fil_WC_8_21_14_2_50_38_16_GCA_002786835.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1                             /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-contigs.fa                             /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG1_02_38_174_GCA_001871945.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1                              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-contigs.fa                              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG2_30_37_12_GCA_001873165.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1            /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-contigs.fa            /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_10_14_0_8_um_filter_38_28_GCA_002785345.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-contigs.fa              /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_CG_4_9_14_0_2_um_filter_38_7_GCA_002788335.1-external-functions.txt
Candidatus_Gracilibacteria_bacterium_GCA_002746815.1                                           /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-contigs.fa                                           /Users/alonshaiber/Downloads/GN02-NCBI-GENOMES/Candidatus_Gracilibacteria_bacterium_GCA_002746815.1-external-functions.txt
```

By excluding the external gene calls file from the fasta.txt file, we guarantee that `anvi-run-workflow` will use prodigal for gene calls.

## Running the workflow

The way I go about this, is that I first run the [contigs workflow](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/#contigs-workflow) on these file.

The easiest thing would be to use this config file `CONTIGS-CONFIG.json`:

```json
{
    "fasta_txt": "GN02-fasta.txt"
}
```

And now we can run:

```bash
anvi-run-workflow -w contigs \
                  -c CONTIGS-CONFIG.json
```

Once this is done, we have generated a contigs database for each of the genomes, and we can add these to an [extrnal genomes file](http://merenlab.org/2016/11/08/pangenomics-v2/#generating-an-anvio-genomes-storage) to use in a pangenomic or phylogenomic analysis.

