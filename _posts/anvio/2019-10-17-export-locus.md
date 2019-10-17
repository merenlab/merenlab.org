# Introduction
Some genetic analyses call for the comparison of specific genetic loci between genomes. For example, one may be interested in investigating evidence for adaptive evolution of the lac operon between different E. coli strains by extracting all loci from different genomes. 

Today, we'll use the tool `anvi-export-locus` to extract the lac operon from a larger genomic context, E. coli genomes!

Briefly, `anvi-export-locus` tool cuts out loci using two approaches: default mode or what we call flank-mode. In the default mode, the tool locates a designated anchor gene, then cuts upstream and downstream context based on user-defined input. Flank-mode, on the other hand, locates designated genes that surround the target locus, then cuts in between them. Target genes of interest to locate anchors for exicion can be defined through their specific ids in anvi'o or through search-terms that query functional annotations or HMM hits stored in your contigs database!

Let's get started.

## Download genomes

Alon has a great tutorial [here](http://merenlab.org/2019/03/14/ncbi-genome-download-magic/) describing how to download genomes from NCBI. We'll be using this to download a few E. Coli genomes today. If you have any questions regarding downloading genomes, please refer to Alon's tutorial.

First, we'll download Genbank files for a few representative E. Coli strains.
```{bash}
# Save working directory for later
WD=$(pwd)

ncbi-genome-download bacteria \
                     --assembly-level chromosome,complete \
                     --genus Escherichia \
                     --metadata metadata.txt \
                     --refseq-category reference
```

This genome is fucked
```{bash}
rm -rf GCF_000008865.2/
```

Nextm we'll make FASTAs, external gene calls, and functional annotations for all the Genbanks we just downloaded.
```{bash}
anvi-script-process-genbank-metadata -m metadata.txt \
                                     --output-dir ecoli \
                                     --output-fasta-txt ecoli.txt
```

# Generate contigs DBs
Now let's run the contigs workflow to creates contigs DBs for each of our genomes/

First, we'll make a json file for the anvio work flow called `contigs.json`
```{bash}
{
    "fasta_txt": "ecoli.txt"
}
```

Then run the workflow!
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

Today with our 5 E. Coli examples we will use a combination of the `default-mode` and `flank-mode` to cut out the lac operon. 

First we will use `default-mode`. This requires you to provide a `-search-term` and `--num-genes` parameter. The `-search-term` will act as an anchor gene to locate the locus within the contigs provided. In this case we will use the lacZ gene (β-galactosidase which cleaves lactose) to locate the lac operon in our  E. coli genomes. Once the anchor gene is located `--num-genes X,Y` will instruct `anvi-export-locus` to cut `X` gene(s) upstream and `Y` gene(s) downstream of the designated anchor gene.

Let's cut of the lac operon to our E. coli genomes using default mode.

The lac operon can be found in different parts of the E. coli genome. 

First we'll use default mode to get the general genomic neighborhood
```{bash}
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

**NOTE**
`--flank-mode` requires flanking genes to be single copies in the contig it's searching. If your locus of interest does not have fixed coordinates in your genomes or metagenomes, you may need to adjust the `-search-terms` on a case by case bases. 

Awesome we have some smaller DNA fragments that contain the lac operon, BUT we also grabbed some extra genes that dont belong to the operon. Let's use `--flank-mode` to trim the loci to just contain the lac operon.

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
