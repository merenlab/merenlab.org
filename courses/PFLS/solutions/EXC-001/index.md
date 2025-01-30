---
layout: page
noleftpanel: true
title: "Solutions for EXC-001"
author: meren
date: February, 2025
geometry: "left=3cm,right=3cm,top=2cm,bottom=2cm"
colorlinks: true
urlcolor: blue
monofont: DejaVuSansMono.ttf
---

[Go back to the questions](../../#exc-001).

## Taking a look at the contents

There are many ways to take a look at the contents of any text file. One can use any of these:

```sh
## see first 20 lines:
head -n 20 genes.fa

## see last 20 lines:
tail -n 20 genes.fa

## browse the contents quickly:
less genes.fa
```

## Finding the number of genes

Even if you didn't know, now you know that each sequence a given FASTA file is marked with a 'defline' that starts with the character `>`. Which means, the number of `>` characters in the file will give us the number of sequences in the FASTA file.


```sh
grep '>' genes.fa | wc -l
```

## Finding the number of genes from a contig

```sh
grep c_000000000023 genes.fa | wc -l
```

## Finding the name of the longest gene

```sh
awk 'BEGIN{FS=":"} />/{print $(NF), $0}' genes.fa | sort -n | tail -n 1 | awk '{print $2}' | sed 's/>//g'
```

## Finding the number genes in a subset of contigs



```sh
awk '/>/{print $2}' genes.fa | awk 'BEGIN{FS=";"} {print $1}' | sort | uniq -c | sort -nr | head -n 10 | awk '{num_genes=num_genes + $1} END {print num_genes}'
```

## Reformatting the FASTA file

So here we wish to turn this FASTA file in which individual sequences can occupy multiple lines and look like this:

```
>0 contig:c_000000000001;start:1580;stop:1754;direction:r;rev_compd:True;length:174
ATGAGAAGTGATGACCCAGTAGGATTAATCGGTGATACCGATAAGAGAACACCATCAGAGATAGCGATTAATGAAGAAGGTAACCATTTGGTAGCTCATAAGCTCGGTGTTATTGAATGC
CCCTATTCTTTAATAGACAAAGCATTAGAGGGCATCTCTAACGTAACCAAATGA
```

into a FASTA file in which each gene sequence fits into a single line like this:

```
>0 contig:c_000000000001;start:1580;stop:1754;direction:r;rev_compd:True;length:174
ATGAGAAGTGATGACCCAGTAGGATTAATCGGTGATACCGATAAGAGAACACCATCAGAGATAGCGATTAATGAAGAAGGTAACCATTTGGTAGCTCATAAGCTCGGTGTTATTGAATGCCCCTATTCTTTAATAGACAAAGCATTAGAGGGCATCTCTAACGTAACCAAATGA
```

There are many ways to do this, but here is how this can happen using this AWK one-liner, which is a great example of 'algorithmic' thinking:

```sh
awk '/>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {print seq}' genes.fa > genes-one-line.fa
```

It is quite difficult, so don't be discouraged if you can't understand it fully just yet.
