---
layout: page
title: A simple read-recruitment exercise
modified: 2020-04-16
excerpt: "With an embarrassingly simple context and reads"
authors: [meren]
comments: true
---

{% include _toc.html %}

The purpose of this read recruitment practice is to get some hands-on experience on technical and rudimentary steps of metagenomic read recruitment and profiling read recruitment results using anvi'o.

## Downloading the data pack

For this exercise we have a little data pack. To download it activate your anvi'o conda environment, open your terminal, and download the package:


```
wget http://merenlab.org/momics/exercises/read-recruitment/momics-week-02-lab-exercise.tar.gz
```

{:.notice}
First run `conda install wget` if you are getting a "command not found" error at this stage. `wget` is a program we commonly use to download things through the terminal. Once `wget` is installed, run the download command again.

Make sure you have the file by running this line:

```
ls -l momics-week-02-lab-exercise.tar.gz
```

If you are not getting an error, you have it :) We first need to unpack it:

```
tar -zxvf momics-week-02-lab-exercise.tar.gz
```

And then go into the new directory that emerged from this:

```
cd momics-week-02-lab-exercise
```

Let's take a look at what we have in the directory:

```
ls

README.md
mock_genome.fa
mock_metagenome_01-R1.fastq
mock_metagenome_01-R2.fastq
mock_metagenome_02-R1.fastq
mock_metagenome_02-R2.fastq
```

Good. We are now ready for the reset of the tutorial.


## Read recruitment practice

This exercise involves two mock metagenomes, called `mock_metagenome_01` and `mock_metagenome_02`, and a genome, `mock_genome`.

What we wish to do is to use the `mock_genome` to recruit reads from our mock metagenomes, and interpret the results.

For read recruitment we will use the program Bowtie2, and if have read [the instructions](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), you know that for mapping to take place, we first need to build a 'reference database' for read recruitment. This is how we can do it for any given reference:

```
bowtie2-build mock_genome.fa mock_genome
```

The next step is to actually align metagenomic short reads against that database (so we would be essentially _recruiting_ reads in these metagenomes that are similar enough to sequences in this database). Let's start with one of our metagenomes, `mock_metagenome_01`:

``` bash
bowtie2 -x mock_genome \
        -1 mock_metagenome_01-R1.fastq \
        -2 mock_metagenome_01-R2.fastq \
        -S mock_metagenome_01.sam
```

The result of this process is a *SAM file*. At this stage you should feel free to take a look at it in a text editor.

What we need, however, is a BAM file, which is the binary form of SAM file. While binary files are not accessible to our human eyes, they are much easier to work with for software since sorting and indexing information in binary files for fast access is much more straightforward.

SAM and BAM are de facto standard file formats in the field as of today, and we can use `samtools` to work with them. Here is how we turn a SAM file into a BAM file, then sort it, and create an index for it:

```
samtools view -F 4 -bS mock_metagenome_01.sam -o mock_metagenome_01-RAW.bam
samtools sort mock_metagenome_01-RAW.bam -o mock_metagenome_01.bam
samtools index mock_metagenome_01.bam
```

The resulting BAM file can be used with many tools. Including [Integrated Genomics Viewer](https://software.broadinstitute.org/software/igv/) or [Artemis BamView](https://www.sanger.ac.uk/science/tools/bamview). They will all show how short reads are stacking underneath the reference context.

But we can also profile them using anvi'o to do even more! But to get anvi'o part of the story going, we first need an *anvi'o contigs database*. An anvi'o contigs database is a special anvi'o file where your reference sequences are characterized.

This is how you generate an anvi'o contigs database:

```
anvi-gen-contigs-database -f mock_genome.fa \
                          -o mock_genome.db
```

Then we can 'profile' our BAM file for the `mock_metagenome_01`:

```
anvi-profile -i mock_metagenome_01.bam \
             -c mock_genome.db \
             -o mock_metagenome_01-profile
```

This command `anvi-profile` creates a **single anvi'o profile database** that enables us to visualize the read recruitment results. Here is how you can view one contig:

```
anvi-inspect -c mock_genome.db \
             -p mock_metagenome_01-profile/PROFILE.db \
             --just-do-it
```

But how about `mock_metagenome_02`? So far we did everything for one of the two metagenomes. To view both of them, we need to do the same steps of read recruitment and profiling for the other metagenome as well. Here are all the steps together:

``` bash
bowtie2 -x mock_genome -1 mock_metagenome_02-R1.fastq -2 mock_metagenome_02-R2.fastq -S mock_metagenome_02.sam
samtools view -F 4 -bS mock_metagenome_02.sam -o mock_metagenome_02-RAW.bam
samtools sort mock_metagenome_02-RAW.bam -o mock_metagenome_02.bam
samtools index mock_metagenome_02.bam
anvi-profile -i mock_metagenome_02.bam -c mock_genome.db -o mock_metagenome_02-profile
```

Now **we have two anvi'o single profile databases**. The next step is to 'merge' these single profiles to have a **merged anvi'o profile database**, so we can see the read recruitment results across multiple samples in the same context. Here is how merging works:

``` bash
anvi-merge mock_metagenome_01-profile/PROFILE.db \
           mock_metagenome_02-profile/PROFILE.db \
           -c mock_genome.db \
           -o merged_profiles
```

If you had many many single profiles to merge, you could have run this command this way:

``` bash
anvi-merge mock_metagenome_*-profile/PROFILE.db \
           -c mock_genome.db \
           -o merged_profiles
```


Now we can take a look at our contig using the merged profile database::

```
anvi-inspect -p merged_profiles/PROFILE.db \
             -c mock_genome.db \
             --just-do-it
```

Tadaa.

{:.warning}
The reason we are using `anvi-inspect` in these examples is because we have only a single tiny contig in our mock genome. In more realistic cases we will use `anvi-interactive` to visualize the distribution of all contigs across all metagenomes. If you want, you can actually replace `anvi-inspect` with `anvi-interactive` to see how it looks like in the last example.


## A more organized way in terminal

The purpose of this section is to mention an aspect of the terminal environment that could make your life easier in the long run.

Let's for a second delete everything we have just generated,

```
rm -rf merged_profiles/ *-profile *bam* *sam* *bt2 *db
```

and talk about the beauties of the shell environment in automatizing repetitive steps of analyses.

First, let's put together all the non-interactive commands we run and stare at them for a second:

``` bash
bowtie2-build mock_genome.fa mock_genome
anvi-gen-contigs-database -f mock_genome.fa -o mock_genome.db
bowtie2 -x mock_genome -1 mock_metagenome_01-R1.fastq -2 mock_metagenome_01-R2.fastq -S mock_metagenome_01.sam
samtools view -F 4 -bS mock_metagenome_01.sam -o mock_metagenome_01-RAW.bam
samtools sort mock_metagenome_01-RAW.bam -o mock_metagenome_01.bam
samtools index mock_metagenome_01.bam
anvi-profile -i mock_metagenome_01.bam -c mock_genome.db -o mock_metagenome_01-profile
bowtie2 -x mock_genome -1 mock_metagenome_02-R1.fastq -2 mock_metagenome_02-R2.fastq -S mock_metagenome_02.sam
samtools view -F 4 -bS mock_metagenome_02.sam -o mock_metagenome_02-RAW.bam
samtools sort mock_metagenome_02-RAW.bam -o mock_metagenome_02.bam
samtools index mock_metagenome_02.bam
anvi-profile -i mock_metagenome_02.bam -c mock_genome.db -o mock_metagenome_02-profile
anvi-merge mock_metagenome_01-profile/PROFILE.db mock_metagenome_02-profile/PROFILE.db -c mock_genome.db -o merged_profiles
```

If you look closely, you can see that certain parts of this workflow include independent steps, while other parts are quite similar and repetitive. Here is an easier to see version of this workflow:

``` bash
# INDEPENDENT / NON-REPETITIVE:
bowtie2-build mock_genome.fa mock_genome
anvi-gen-contigs-database -f mock_genome.fa -o mock_genome.db

# REPETITIVE PART 1:
bowtie2 -x mock_genome -1 mock_metagenome_01-R1.fastq -2 mock_metagenome_01-R2.fastq -S mock_metagenome_01.sam
samtools view -F 4 -bS mock_metagenome_01.sam -o mock_metagenome_01-RAW.bam
samtools sort mock_metagenome_01-RAW.bam -o mock_metagenome_01.bam
samtools index mock_metagenome_01.bam
anvi-profile -i mock_metagenome_01.bam -c mock_genome.db -o mock_metagenome_01-profile

# REPETITIVE PART 2:
bowtie2 -x mock_genome -1 mock_metagenome_02-R1.fastq -2 mock_metagenome_02-R2.fastq -S mock_metagenome_02.sam
samtools view -F 4 -bS mock_metagenome_02.sam -o mock_metagenome_02-RAW.bam
samtools sort mock_metagenome_02-RAW.bam -o mock_metagenome_02.bam
samtools index mock_metagenome_02.bam
anvi-profile -i mock_metagenome_02.bam -c mock_genome.db -o mock_metagenome_02-profile

# INDEPENDENT / NON-REPETITIVE:
anvi-merge mock_metagenome_01-profile/PROFILE.db mock_metagenome_02-profile/PROFILE.db -c mock_gen=ome.db -o merged_profiles
```

The repetitive parts actually follow this template:

``` bash
bowtie2 -x mock_genome -1 ${metagenome}-R1.fastq -2 ${metagenome}-R2.fastq -S ${metagenome}.sam
samtools view -F 4 -bS ${metagenome}.sam -o ${metagenome}-RAW.bam
samtools sort ${metagenome}-RAW.bam -o ${metagenome}.bam
samtools index ${metagenome}.bam
anvi-profile -i ${metagenome}.bam -c mock_genome.db -o ${metagenome}-profile
```

So, this template would have run everything successfully in the *REPETITIVE PART 1* if were to replace all `${metagenome}` values with `mock_metagenome_01`. Similarly, it would work also for *REPETITIVE PART 2* if were to replace all `${metagenome}` values with `mock_metagenome_02`.

That is exactly shell scripting enables us to do. The following BASH script will do exactly what we have done before with much less number of lines (and less room for human error):


``` bash
# INDEPENDENT / NON-REPETITIVE:
bowtie2-build mock_genome.fa mock_genome
anvi-gen-contigs-database -f mock_genome.fa -o mock_genome.db


# REPETITIVE PART for two samples:
for metagenome in mock_metagenome_01 mock_metagenome_02
do
    bowtie2 -x mock_genome -1 ${metagenome}-R1.fastq -2 ${metagenome}-R2.fastq -S ${metagenome}.sam
    samtools view -F 4 -bS ${metagenome}.sam -o ${metagenome}-RAW.bam
    samtools sort ${metagenome}-RAW.bam -o ${metagenome}.bam
    samtools index ${metagenome}.bam
    anvi-profile -i ${metagenome}.bam -c mock_genome.db -o ${metagenome}-profile
done

# INDEPENDENT / NON-REPETITIVE:
anvi-merge mock_metagenome_01-profile/PROFILE.db mock_metagenome_02-profile/PROFILE.db -c mock_genome.db -o merged_profiles
```

Now you can store this into a file, let's call it `exercise.sh`, and run it like this to get all these steps done:

```
bash exercise.sh
```

and visualize the results when everything is finished:

```
anvi-inspect -p merged_profiles/PROFILE.db \
             -c mock_genome.db \
             --just-do-it
```

It looks like a small improvement in this example. But it would have mattered much more if you had many more samples to process :)

**The next stage of automatizing and managing large number of jobs** is to use workflow building software, such as Snakemake or NextFlow. Anvi'o has an infrastructure for that which enables us to work with thousands of samples and billions of metagenomic reads without even thinking about them. If you are interested in taking a look, you can learn about anvi'o snakemake workflows [here](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/).