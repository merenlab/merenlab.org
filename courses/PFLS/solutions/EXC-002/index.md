---
layout: page
noleftpanel: true
title: "Solution for EXC-002"
author: meren
date: February, 2025
geometry: "left=3cm,right=3cm,top=2cm,bottom=2cm"
colorlinks: true
urlcolor: blue
monofont: DejaVuSansMono.ttf
---

[Go back to the question](../../#exc-002).

```sh
## you haven't learned about this, but this is a special instruction for a BASH script
## to stop when there is an error, instead of going on and on
set -e

## set the maximum number of contigs desired
max_num_contigs=$1

## define a variable for the output directory
output_dir="WOLBACHIA-GENOMES"

## create the output directory, if it is not there.
if [ ! -d $output_dir ]; then
    mkdir $output_dir
fi

## go through each genome, and do everything all at once:
for fasta_file in *fna; do
    # learn the number of contigs in the FASTA file
    num_contigs=$(grep '>' $fasta_file | wc -l)

    if [ $num_contigs -lt $max_num_contigs ]; then
        # if we are here, it means the genome has less number of
        # contigs than the threshold defined by the user. now it
        # is time to determine the genome name from the FASTA
        # file
        genome_name=$(echo $fasta_file | awk 'BEGIN {FS="_"} {print $1 "_" $2}')

        # grep the information from the wolbachia-hosts.txt using $genome_name
        new_fasta_name=$(grep $genome_name wolbachia-hosts.txt | awk '{print $2}')

        # copy the $fasta_file for this genome into the $output_dir with its $new_fasta_name
        cp $fasta_file $output_dir/$new_fasta_name.fa
    fi
done
```
{% include CODEBLOCKFILENAME filename="exc-002-eren.sh" %}
