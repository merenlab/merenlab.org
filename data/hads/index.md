---
layout: page
title: Hawaiʻi Diel Sampling (HaDS)
modified: 2024-12-18
excerpt: ""
comments: true
authors: []
---

The purpose of this page is to provide access to the raw data and reproducible data products generated from the Hawaiʻi Diel Sampling (HaDS) Project.

We are still in the process of preparing and updating the contents of this page. Please keep an eye on this space for more soon.

{% include IMAGE path="images/hawaii-diel-2021-74.jpg" width=80 caption="Transiting to our offshore sampling site through the Sampan Channel of Kāneʻohe Bay, Oʻahu, Hawaiʻi on the first morning of sample collection." %}

## Motivation

Microbial communities experience environmental fluctuations across timescales ranging from seconds to seasons, and their responses are evident at multiple levels -- from changes in community composition to the physiological reactions of individual cells or from diel cycles to seasonal variations ([Fuhrman et al. 2015](https://www.nature.com/articles/nrmicro3417)).

Time-series studies in marine systems have largely focused on resolving changes in microbial community composition from seasonal ([Bunse and Pinhassi 2017](https://doi.org/10.1016/j.tim.2016.12.013); [Giovannoni and Vergin 2012](https://doi.org/10.1126/science.1198078)) to daily timescales ([Needham and Fuhrman 2016](https://doi.org/10.1038/nmicrobiol.2016.5)). While the physiological responses of individual microbial populations offer important insights into their ecology and evolution, such population level responses, especially at short timescales, are less well understood in complex environments. Responses to short-term fluctuations that occur on timescales that span from seconds to hours are mostly reflected in changes at the level of transcription and translational regulation without any immediate impact on community composition. We generated the HaDS dataset to contribute an interlinked 'omics resource that lends itself to studies of subtle and population-resolved responses of microbes to environmental variability.

HaDS is a collection of metagenomes, metatranscriptomes and metaepitranscriptomes generated over a 48 hour period at 90 minutes intervals at two sampling sites in Kāneʻohe Bay, Hawaiʻi. The spatiotemporal dynamics of the two surface water sampling stations (so-called HP1 and STO1) are well characterized through the Kāneʻohe Bay Time-series ([Tucker et al. 2021](https://doi.org/10.7717/peerj.12274)), an ongoing monthly time-series sampling program of surface ocean biogeochemistry and microbial communities. Our high-resolution multi-omics approach, paired with concurrent measurements of biogeochemical parameters (chlorophyll, temperature, and nutrient concentration) and contextualized by long-term microbial community and biogeochemistry data at both sampling sites, enables the exploration of microbial population responses to environmental fluctuations and long-term change.

{% include IMAGE path="images/hawaii-diel-2021-80.jpg" width=80 caption="Sample processing next to the docks of Hawaiʻi Institute of Marine Biology (HIMB)." %}

## Data

At both the coastal Kāneʻohe Bay station (HP1) and the adjacent offshore station (STO1), we sampled at 33 time-points across 48 hours. We subsequently produced 59 metatranscriptomes, 65 short-read metagenomes, 8 long-read metagenomes, and 66 metaepitranscriptomes. We also generated four deeply-sequenced short-read metagenomes from samples collected in the late fall and spring prior to HaDS through routine Kāneʻohe Bay Time-series sampling. The following data items give access to RAW sequencing results  as well as processed data items through repositories at [FigShare](https://figshare.com/projects/Hawai_i_Diel_Sampling_HADS_/244907), BCO-DMO, and NCBI:

* NCBI Project ID [PRJNA1201851](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1201851) offers access to [all raw data for short-read and long-read metagenomes, as well as metatranscriptomes and metaepitranscriptomes](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1201851).
* doi:(pending URL from BCO-DMO) provides access to all biogeochemical data that covers the sampling period.
* doi:[10.6084/m9.figshare.28784717](https://doi.org/10.6084/m9.figshare.28784717) serves anvi'o {% include ARTIFACT name="contigs-db" %} files for the individual co-assemblies of short-read (SR) as well as long-read (LR) sequencing of metagenomes. Please note that an anvi'o {% include ARTIFACT name="contigs-db" %} includes gene calls, functional annotations, HMM hits, and other information about each contig, and you can always use the program {% include PROGRAM name="anvi-export-contigs" %} to get a FASTA file for sequences. The following screenshot of the {% include PROGRAM name="anvi-display-contigs-stats" %} output gives an idea about the contents of each {% include ARTIFACT name="contigs-db" %} file:

{% include IMAGE path="images/anvi-display-contigs-stats.png" width=80 caption="A screenshot of the `anvi-display-contigs-stats` output." %}

* doi:[10.6084/m9.figshare.28784762](https://doi.org/10.6084/m9.figshare.28784762) serves FASTA files for metagenome-assembled genomes (MAGs) we have reconstructed from short-read and long-read sequencing of the metagenomes. They are the outputs of quite a preliminary effort, thus secondary attempts to recover genomes from the co-assemblies are most welcome (and very much encouraged). Please see the Supplementary Table for taxonomic annotation and completion / redundancy estimates of the MAGs.
* doi:[10.6084/m9.figshare.28784765](https://doi.org/10.6084/m9.figshare.28784765). The [EcoPhylo](https://anvio.org/help/main/workflows/ecophylo/) output that describes the phylogeography of ribosomal protein L14.

Please feel free to reach out to us if you have any questions regarding access and/or processing of these datasets.

## Bioinformatics

The purpose of this section is to describe key steps of our data generation and formatting workflow, products of which are shared in the previous section. For all downstream analyses, we used long-read metagenomes, quality-filtered short-read metagenomes, and adapter-trimmed and quality-filtered short-read metatranscriptomes. tRNA sequencing data required custom steps of demultiplexing which we describe in greater detail later in this document.

{:.notice}
The primary purpose of the following commands is to give the reader an overall understanding of the bioinformatics steps rather than offering a truly reproducible recipe. In our analyses, we used our high-performance computing clusters to parallelize many of these steps. Please feel free to reach out to us if you have any questions.

We then co-assembled short-read and long-read sequencing data from each station using IDBA-UD and hifiasm-meta, respectively.

### Processing of co-assemblies, read mapping, and binning

We generated the anvi'o {% include ARTIFACT name="contigs-db" %} files that are stored at doi:[10.6084/m9.figshare.28784717](https://doi.org/10.6084/m9.figshare.28784717) using the following commands:

```bash
num_threads="40"
for station in xHP1 STO1
do
    for technology in SR LR
    do
        # generate contigs-db file
        anvi-gen-contigs-database -f ${station}-${technology}-COASSEMBLY.fa -o ${station}-${technology}-COASSEMBLY.db -T ${num_threads}

        # annotate genes
        anvi-run-hmms -c ${station}-${technology}-COASSEMBLY.db -T ${num_threads}
        anvi-scan-trnas -c ${station}-${technology}-COASSEMBLY.db -T ${num_threads}
        anvi-run-ncbi-cogs -c ${station}-${technology}-COASSEMBLY.db -T ${num_threads}
        anvi-run-scg-taxonomy -c ${station}-${technology}-COASSEMBLY.db -T ${num_threads}
    done
done
```

We then used Bowtie2 the following way to recruit short metagenomic and metatranscripomic reads, which were stored in a {% include ARTIFACT name="samples-txt" %} file called `samples-txt.txt` in our working directory, to all assemblies separately and profiled them using anvi'o:


```bash
num_threads="40"
for station in xHP1 STO1
do
    for technology in SR LR
    do
        # generate a directory to store mapping results for each co-assembly
        mkdir ${station}-${technology}

        # generate a bowtie index
        bowtie2-build ${station}-${technology}-COASSEMBLY.fa ${station}-${technology}/${station}-${technology}-COASSEMBLY

        while read sample r1 r2;
        do
            if [ "$sample" == "sample" ]; then continue; fi

            # generate the SAM file
            bowtie2 --threads ${num_threads} \
                    -x ${station}-${technology}/${station}-${technology}-COASSEMBLY \
                    -1 $r1 \
                    -2 $r2 \
                    --no-unal \
                    -S ${station}-${technology}/$sample.sam

            # covert the resulting SAM file to a BAM file:
            samtools view -F 4 -bS ${station}-${technology}/${sample}.sam > ${station}-${technology}/${sample}-RAW.bam

            # sort and index the BAM file:
            samtools sort ${station}-${technology}/$sample-RAW.bam -o ${station}-${technology}/$sample.bam
            samtools index ${station}-${technology}/$sample.bam

            # remove temporary files:
            rm ${station}-${technology}/$sample.sam ${station}-${technology}/$sample-RAW.bam
        done < samples-txt.txt

        # mapping for ${station}-${technology} is done, now we can profile the resulting
        # BAM files
        while read sample r1 r2;
        do
            if [ "$sample" == "sample" ]; then continue; fi

            anvi-profile -c ${station}-${technology}-COASSEMBLY.db \
                         -i ${station}-${technology}/$sample.bam \
                         --profile-SCVs \
                         -M 100 \
                         --num-threads ${num_threads} \
                         -o ${station}-${technology}/$sample
        done < samples-txt.txt

        # all single profiles are ready, and now we can merge them to get a
        # merged anvio profile
        anvi-merge ${station}-${technology}/*/PROFILE.db -o ${station}-${technology}-MERGED -c ${station}-${technology}-COASSEMBLY.db
    done
done
```

We used metabat2 to bin contigs and imported the resulting bins into anvi'o the anvi'o merged {% include ARTIFACT name="profile-db" %} as a {% include ARTIFACT name="collection" %}:

```bash
conda activate metabat2

mkdir HADS-MAGs

for station in xHP1 STO1
do
    # identify BAM files that belong to the diel sampling and describe short-read
    # metagenomic read recruitment results
    BAM_FILES=$(ls ${station}-SR/HADS_202108*MGX*.bam)

    jgi_summarize_bam_contig_depths --outputDepth ${station}_depth.txt --pairedContigs ${station}_paired.txt $BAM_FILES

    metabat2 -t ${num_threads} -i ${station}-SR-COASSEMBLY.fa -a ${station}_depth.txt -o HADS-MAGs/${station}-SR-COASSEMBLY-BIN -v

    FILES=$(ls HADS-MAGs/{station}-SR-COASSEMBLY-BIN*.fa)
    for f in $FILES
    do
        NAME=$(basename $f .fa)
        grep ">" $f | sed 's/>//' | sed -e "s/$/\t$NAME/" | sed 's/\./_/' >> ${station}-collection.txt
    done

    # finally, import the collection into the merged profile:
    anvi-import-collection ${station}-collection.txt -p ${station}-SR-MERGED/PROFILE.db -c ${station}-SR-COASSEMBLY.db
done
```

### EcoPhylo analysis

We followed the [EcoPhylo](https://anvio.org/help/main/workflows/ecophylo/) workflow to characterize the biogeography of ribosomal protein L14 sequences in our data, and to manually curate the ribosomal protein phylogeny by removing sequences that were not classified at the domain level and found in branches that primarily described chloroplast or mitochondrial genomes.

### Demultiplexing of tRNA sequencing results

The tRNA sequencing uses specific barcodes to multiplex samples during library preparation. We used the following script to demultiplex the raw sequencing data prior to uploading them to the NCBI:

```python
# -*- coding: utf-8 -*-

import argparse
import gzip
import os


def open_fastq(file_path):
    """Open a FASTQ file, using gzip if necessary."""
    return gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path)


def create_output_files(output_dir, samples, barcodes):
    """Create output FASTQ files for each barcode-sample pair."""
    r1_outputs = {}
    r2_outputs = {}
    for sample, barcode in zip(samples, barcodes):
        r1_path = os.path.join(output_dir, f"{sample}.r1.fastq.gz")
        r2_path = os.path.join(output_dir, f"{sample}.r2.fastq.gz")
        r1_outputs[barcode] = gzip.open(r1_path, 'wt')
        r2_outputs[barcode] = gzip.open(r2_path, 'wt')
    return r1_outputs, r2_outputs


def parse_args():
    parser = argparse.ArgumentParser(description='Demultiplex paired-end FASTQ files by barcode.')
    parser.add_argument('--r1', required=True, help='Read 1 FASTQ file')
    parser.add_argument('--r2', required=True, help='Read 2 FASTQ file')
    parser.add_argument('--location', choices=('r1', 'r2'), required=True,
                        help='Location of barcode: beginning of Read 1 or Read 2')
    parser.add_argument('--barcodes', nargs='+', required=True, help='List of barcode sequences')
    parser.add_argument('--samples', nargs='+', required=True, help='Sample names corresponding to barcodes')
    parser.add_argument('--outdir', required=True, help='Directory for output files')
    return parser.parse_args()


def demultiplex(r1_file, r2_file, barcodes, r1_outputs, r2_outputs, barcode_location):
    """Demultiplex the FASTQ read pairs based on barcodes."""
    barcode_in_r1 = (barcode_location == 'r1')
    seq_line_counter = 0
    r1_block = []
    r2_block = []
    matched_barcode = None

    for r1_line, r2_line in zip(r1_file, r2_file):
        seq_line_counter += 1
        r1_block.append(r1_line)
        r2_block.append(r2_line)

        if seq_line_counter == 2:
            sequence_line = r1_line if barcode_in_r1 else r2_line
            for barcode in barcodes:
                if sequence_line.startswith(barcode):
                    matched_barcode = barcode
                    break

        if seq_line_counter == 4:
            if matched_barcode:
                r1_outputs[matched_barcode].writelines(r1_block)
                r2_outputs[matched_barcode].writelines(r2_block)
            # Reset for next block
            r1_block = []
            r2_block = []
            seq_line_counter = 0
            matched_barcode = None


def main():
    args = parse_args()

    with open_fastq(args.r1) as r1_file, open_fastq(args.r2) as r2_file:
        r1_outputs, r2_outputs = create_output_files(args.outdir, args.samples, args.barcodes)
        try:
            demultiplex(r1_file, r2_file, args.barcodes, r1_outputs, r2_outputs, args.location)
        finally:
            for output_file in list(r1_outputs.values()) + list(r2_outputs.values()):
                output_file.close()


if __name__ == '__main__':
    main()
```
