---
layout: post
title: "Importing VirSorter annotations into anvi'o to study phages"
excerpt: "Anvi'o projects meet the underappreciated owners of this planet."
modified: 2018-02-08
tags: [tutorial, Virsorter]
categories: [anvio]
comments: true
authors: [bryan]
---

{% capture images %}{{site.url}}/images/anvio/2018-02-08-working-with-VirSorter{% endcapture %}

{% include _toc.html %}

{:.notice}
**A note from the Meren Lab**: We are very thankful to Bryan for sharing his expertise on behalf of all anvi’o users who wished to know more about the phages lurking in their metagenomes using [VirSorter](https://doi.org/10.7717/peerj.985). Bryan is currently a PhD Candidate in the [Sonnenburg Lab](http://sonnenburglab.stanford.edu/), Stanford University School of Medicine.

Metagenome studies are a great way to explore complex communities. Many algorithms and tools have made it possible to reconstruct bacterial, archaeal, and even eukaryotic genome bins for diverse organisms from assembled metagenome sequencing data. However, these algorithms are not well equipped to deal with the most abundant biological entity on earth -- bacteriophages. There are several tools designed to predict viral contigs from metagenome assemblies, and prophages from bacterial genomes.

This tutorial will walk you through the steps needed to (1) use VirSorter to predict which contigs in your assembly are phages, and (2) visualize these results in [anvi'o](https://doi.org/10.7717/peerj.1319) using the programs `anvi-interactive` and `anvi-refine`. The only thing better than binning with anvi'o is _phage-aware_ binning with anvi'o. Here we go!


## A brief introduction to VirSorter

VirSorter was published in [Roux et al (2015)](https://doi.org/10.7717/peerj.985 "VirSorter: mining viral signal from microbial genomic data"), along with its [source code](https://github.com/simroux/VirSorter). The following figure ([Figure 1A](https://doi.org/10.7717/peerj.985/fig-1) in the original study) explains the VirSorter pipeline:

[![virsorter](https://dfzljdn9uc3pi.cloudfront.net/2015/985/1/fig-1-2x.jpg)](https://dfzljdn9uc3pi.cloudfront.net/2015/985/1/fig-1-2x.jpg){:.center-img .width-50}

Briefly:

- The VirSorter input is a single FASTA file. In our case this will be the same metagenome assembly you used to make your anvi'o contigs database.

- VirSorter annotates the FASTA file using MetaGeneAnnotator ([Noguchi et al, 2006](https://doi.org/10.1093/nar/gkl723)), and then uses hmmsearch ([Eddy et al, 2011](https://doi.org/10.1371/journal.pcbi.1002195)) to predict (1) PFAM domains (Version 27, [Finn et al, 2014](https://doi.org/10.1093/nar/gkt1223)) and (2) viral domains using pre-computed, downloadable HMM databases on the annotated genes.

- For each contig, VirSorter uses a sliding window analysis to identify regions of several genes that: (1) contain one or more viral "hallmark" genes (capsid, large terminase, etc.), (2) are enriched for viral domains, (3) have few PFAM domains, (4) have many uncharacterized genes, (5) have many short genes, and (6) have many genes encoded sequentially on the same strand. Overlapping gene regions predicted by each criterion are combined. Depending on which criteria are met, the contig is assigned a category number, where 1 is high confidence that it the contig is a phage, and 3 is possible, but low confidence.

- If a category 1, 2, or 3 prediction on a conting encompasses > 80% of the contig, the whole contig is annotated as being a "phage". If a category 1, 2, or 3 prediction on a contig encompasses <= 80% of the contig, then a subset of the contig is annotated as being a "prophage". For metagenome assemblies, note that a contig annotated as being a "phage" could actually be part of a prophage but isn't called as such simply becasue more than 80% of the genes are "phage-like".

## Installing VirSorter locally with conda

If you don't want to install anything, you can run VirSorter on [CyVerse](http://user.cyverse.org/). Alternatively, you can run VirSorter locally [using Docker](https://github.com/simroux/VirSorter), or by installing the codebase on your computer. Here I will explain how to install the VirSorter codebase on your own computer, which is what I recommend you consider doing.

To manually install the VirSorter codebase you can follow the recipe down below. I have only tested these instructions on Linux (Ubuntu and CentOS), but if you have any questions feel free to contact me or leave a comment down below.

Assuming you have [conda](https://conda.io/docs/user-guide/install/index.html) installed, copy-paste the following commands in your terminal after going into the directory where you wish to install VirSorter:

``` bash
# create a conda environment
conda create --name virsorter \
             -c bioconda mcl=14.137 \
                         muscle \
                         blast \
                         perl-bioperl \
                         perl-file-which \
                         hmmer=3.1b2 \
                         perl-parallel-forkmanager \
                         perl-list-moreutils

# clone the VirSorter repository and go into it
git clone https://github.com/simroux/VirSorter.git
cd VirSorter
git checkout e98d2f8f473b3028793b8a037c91648d1453f7a0 #Checks out version 1.0.5
cd Scripts

# run make
make
```

{:.notice}
VirSorter v.1.0.5 is much faster than previous versions because several steps are now parallelized. Unless you opt to use diamond (see VirSorter usage), the results from running VirSorter v1.0.5 should be identical to running v1.0.3.

The commands above will create everything you need to run VirSorter, but VirSorter commands will not be available to you systemwide. To run VirSorter from any directory, you can make symbolic links to `VirSorter/wrapper_phage_contigs_sorter_iPlant.pl` and `VirSorter/Scripts` and place them in the `bin` folder for your "virsorter" conda environment. An example location of this `bin` folder is `~/miniconda/envs/virsorter/bin`. Substitute this path with the path to the `bin` folder for your newly created "virsorter" environment.

``` bash
ln -s ~/Applications/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl ~/miniconda/envs/virsorter/bin

ln -s ~/Applications/VirSorter/Scripts ~/miniconda/envs/virsorter/bin
```

To run VirSorter you will also need to download MetaGeneAnnotator ([Noguchi et al, 2006](https://doi.org/10.1093/nar/gkl723)). I like to just put this in the virsorter environment's `bin` folder alongside the VirSorter symbolic links.

```
cd ~/miniconda/envs/virsorter/bin
wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
tar -xvzf metagene/mga_x86_64.tar.gz
```

Now you have a running VirSorter installation, however, you need the data pack to actually run it.


## Obtaining the VirSorter data pack

The virsorter-data archive contains:

1. PFAM-A and PFAM-B HMM models (Version 27, [Finn et al, 2014](https://doi.org/10.1093/nar/gkt1223))

2. Two HMM databases computed from (1) all phages in RefSeq prior to 2014, and (2) those same phages plus curated phages from several viromes

3. Other files VirSorter needs in order to run

Navigate to a directory where you want the data pack to live, and run the following commands. It doesn't have to be the same location where you downloaded VirSorter.

```
# download the data pack I re-formatted from the original data pack
wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz

# run md5sum to make sure the output
# matches to dd12af7d13da0a85df0a9106e9346b45
md5sum virsorter-data-v2.tar.gz

# unpack it
tar -xvzf virsorter-data-v2.tar.gz
```

To run VirSorter, you will need to point to the path of the data-pack directory using the `--data-dir` parameter each time you run VirSorter.

<div class="extra-info" markdown="1">

<span class="extra-info-header">Why a re-formatted data-pack?</span>

The original HMMs made available by Roux et al (2015) on CyVerse were formatted for HMMER version 3.0 (and were very hard to download). The HMMs in the updated data-pack are compatible with HMMER version 3.1b2.

In order to generate this updated datapack, I first downloaded the original data-pack from CyVerse, and then run the following command in it to re-format it:

``` bash
for i in */*.hmm
do
    echo "Converting ${i}..."
    hmmconvert ${i} > ${i}.new
    mv ${i}.new ${i}
    hmmpress -f ${i}
done
```

Then Meren uploaded this datapack [here](https://doi.org/10.5281/zenodo.1168727) so it can be accessed publicly.

</div>


## Running VirSorter

You can always run VirSorter on CyVerse or using Docker. But why not just install it yourself? Once you have the conda environment created, VirSorter downloaded, and the symbolic links made, running VirSorter is as easy as:

```
source activate virsorter
wrapper_phage_contigs_sorter_iPlant.pl -f assembly.fasta \
                                       --db 1 \
                                       --wdir output_directory \
                                       --ncpu 4 \
                                       --data-dir /path/to/virsorter-data
```

{:.warning}
It is _critical_ that VirSorter is run on the _exact same FASTA file_ that was used to generate the anvi'o contigs database, map and profile reads, etc. This point can't be overstated!

Please note that the `--db` argument can either be `1` or `2`. If set to `1`, then VirSorter will use phage HMMs computed from RefSeq phages published before January 2014. If set to `2`, then Virsorter will use all of the HMMs from option `1`, plus additional HMMs from phage contigs identified in curated virome datasets.


## Gathering files needed to run VirSorterParser

VirSorter writes several outputs to a working directory that you specify when run VirSorter. The files we need for importing VirSorter annotations into anvi'o include:

- `VIRSorter_global-phage-signal.csv`: This file contains one line for each phage prediction. Often, this results in one line per contig, though very large contigs might have two or more prophage predictions.

- `Metric_files/VIRSorter_affi-contigs.tab`: This file contains one line per gene and includes any PFAM or phage domain annotations. The lines for all genes for a given contig are preceded by a line containing `>Contig_name`.

If you have been using anvi'o for your metagenomes, I assume you already have your anvi'o contigs database.

Run the following command to generate one of the essential files we will need to visualize VirSorter split-level annotations in anvi'o:

``` bash
anvi-export-table CONTIGS.db \
                  --table splits_basic_info
```

This will generate an output file called `splits_basic_info.txt` file, which we will use alongside with the output files generated VirSorter to import everything into anvi'o.

If you want to also generate gene-level annotations for hallmark genes and import those into anvi'o (think capsid proteins, terminases, etc.) then you'll need to also run this command to export your gene calls from your anvi'o contigs database:

``` bash
anvi-anvi-export-gene-calls CONTIGS.db \
                  -o all_gene_calls.txt
```


To make sense of all these files, I implemented a program called `VirSorterParser`.

<div class="extra-info" markdown="1">

<span class="extra-info-header">Rationale for the VirSorterParser</span>

VirSorter uses MetaGeneAnnotator ([Noguchi et al, 2006](https://doi.org/10.1093/nar/gkl723)) and doesn't accept custom gene calls. Some of the predicted genes may not line up with what anvi'o has predicted on your contigs database. However, if you export your anvi'o gene calls and pass that file to the parser, the parser will attempt to get around this by matching up the start/stop positions of hallmark genes predicted by VirSorter with existing gene calls from your anvi'o contigs database. For each hallmark gene in a contig identified as phage or prophage where this matching is successful, the parser will write that gene and the VirSorter-predicted function to a file that can be imported into your contigs database using `anvi-import-functions`. 

The parser generates an additional data file that can be visualized when running anvi-interactive or anvi-refine. The additional data file the parser will generate follows this structure:

|split|phage_name|category|num_genes_in_phage|num_phage_hallmark_genes_in_phage|
|:--|:--|:--|:--|:--|
|(...)|(...)|(...)|(...)|(...)|

For each VirSorter phage or prophage prediction that spans several splits, all other columns besides "split" are identical across splits. These four metrics are reported on the predicted phage or prophage, not on a given split.

For example, if there are 86 genes and 3 splits in the phage contig, `num_genes_in_phage` will report `86` for `split_00001`, `split_00002`, and `split_00003`. The same is true for `num_phage_hallmark_genes_in_phage`.

Reported categories in the `category` column include the following (see [Figure 1B](https://doi.org/10.7717/peerj.985/fig-1) original publication, and information therein):

- cat1_phage
- cat2_phage
- cat3_phage
- cat1_prophage
- cat2_prophage
- cat3_prophage

For the `phage_name` column, the first phage predicted by VirSorter is named `phage_1` and increments by one up to `phage_n`. Similarly, the first prophage predicted by VirSorter is named `prophage_1` and incrememts similarly.
</div>

### Running the parser

To run the parser you just need the python script `virsorter_to_anvio.py`. If you want to import hallmarn gene functions into anvi'o, you'll need two mapping files (one for each VirSorter database). You can download these files into your work directory the following way:

``` bash
wget https://raw.githubusercontent.com/brymerr921/VirSorterParser/master/virsorter_to_anvio.py
wget https://raw.githubusercontent.com/brymerr921/VirSorterParser/master/hallmark_to_function_files/db1_hallmark_functions.txt
wget https://raw.githubusercontent.com/brymerr921/VirSorterParser/master/hallmark_to_function_files/db2_hallmark_functions.txt
```

Arguments are described in the help menu:

``` bash
python virsorter_to_anvio.py --help
```

Briefly, the script takes as input the two output files from VirSorter, `VIRSorter_affi-contigs.tab` and `VIRSorter_global_signal.csv`, and the `splits_basic_info.txt` file from anvi'o. It also asks you which VirSorter database you used when you ran VirSorter. These are required inputs. If you want to test out the parser, samples of each of the required files are provided in the "sample_files" directory of [this repository](https://github.com/brymerr921/VirSorterParser). If you run parser on the the sample files, specify "--db 2" on the command line.

You can control which VirSorter predictions are prepared for importing into anvi'o:

* `--exclude-cat3` will skip over all "category 3" predictions.

* `--exclude-prophages` will skip over all prophages.

* The `-L` parameter can be used to specify the minimum phage length to report in the output files. For example, `-L 5000` means that all phage predictions shorter than 5000 bp will be not reported in the output files. These flags can be used in any combination with each other.

If you want to import hallmark gene functions from predicted phage or prophage contigs into anvi'o, you'll need to also supply the parser with `--anvio-gene-calls all_gene_calls.txt` and also a mapping file available in the [VirSorterParser repository](https://github.com/brymerr921/VirSorterParser/tree/master/hallmark_to_function_files): `--hallmark-functions hallmark_to_functions_files/db2_hallmark_functions.txt`.

An example command to run the parser is:

``` bash
./virsorter_to_anvio.py -a VIRSorter_affi-contigs.tab -g VIRSorter_global-phage-signal.csv -s splits_basic_info.txt --db 2 --anvio-gene-calls all_gene_calls.txt --hallmark-functions db2_hallmark_functions.txt --db 2
```

### Output files

The parser will generate two output files by default. The first is "virsorter_additional_info.txt", which can be used direcly as an additional data file with anvi'o.

For anvi'o `v3` and earlier, you use this file when running `anvi-interactive` or `anvi-refine` by using the flag `-A virsorter_additional_info.txt`.

For anvi'o `v4` or later, you can import this file into your profile database as [additional data](http://merenlab.org/2017/12/11/additional-data-tables/):

``` bash
anvi-import-misc-data virsorter_additional_info.txt \
                      -p PROFILE.db \
                      --target-data-table items
```

The second output file is "virsorter_collection.txt". This file can be imported into anvi'o as a collection:

``` bash
anvi-import-collection virsorter_collection.txt\
                       -c CONTIGS.db \
                       -p PROFILE.db \
                       -C COLLECTION_NAME
```

This will generate a collection where each phage prediction will become a bin containing all of the splits for that phage.


If you specified files for `--anvio-gene-calls` and `--hallmark-functions`, you'll also get an output file called "virsorter_annotations.txt" containing gene annotations for all hallmark genes inside all contigs predicted as phage or prophage by VirSorter. You can import it into anvi'o like this:

``` bash
anvi-import-functions 
                       -c CONTIGS.db \
                       -i virsorter_annotations.txt
```

## Happy binning!

Congratulations, you can now enjoy an even better _phage-aware_ binning experience in anvi'o! If you have any questions or problems, please don't hesitate to contact me, or submit an issue on my [GitHub repository](https://github.com/brymerr921/VirSorterParser).

