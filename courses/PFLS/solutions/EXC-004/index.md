---
layout: page
noleftpanel: true
title: "Solution for EXC-004"
author: meren
date: February, 2025
geometry: "left=3cm,right=3cm,top=2cm,bottom=2cm"
colorlinks: true
urlcolor: blue
monofont: DejaVuSansMono.ttf
---

[Go back to the question](../../#exc-004).

```sh
## remove the COMBINED-DATA directory if it exists from a previous run
rm -rf COMBINED-DATA

## create a fresh one -- by now we know that there is one and
## only COMBINED-DATA directory and it is empty
mkdir -p COMBINED-DATA

## go through all directories that have FASTA files for different
## cultures
for dir in $(ls -d RAW-DATA/DNA*); do
    # learn the culture name from the directory path:
    culture_name=$(basename $dir)

    # learn the new name for a given culture from the
    # sample-translation.txt file
    new_culture_name=$(grep $culture_name RAW-DATA/sample-translation.txt | awk '{print $2}')

    # since we want to number final FASTA files nicely after figuring out
    # what is a bin and what is a MAG based on the completion estimates,
    # we need variables that keep track of those numbers
    MAG_counter=1
    BIN_counter=1

    # get the files that include completion estimates and taxonomy
    # into the new output directory
    cp $dir/checkm.txt COMBINED-DATA/$new_culture_name-CHECKM.txt
    cp $dir/gtdb.gtdbtk.tax COMBINED-DATA/$new_culture_name-GTDB-TAX.txt

    # this is our SECOND for loop to go through each FASTA file in the
    # input directory
    for fasta_file in $dir/bins/*.fasta; do
        # determine the bin name from the file name itself
        bin_name=$(basename $fasta_file .fasta)

        # learn the completion and contamination estimates for the bin
        completion=$(grep "$bin_name " $dir/checkm.txt | awk '{print $13}')
        contamination=$(grep "$bin_name " $dir/checkm.txt | awk '{print $14}')

        # determine if we should name it with the postfix UNBINNED, BIN, or MAG, and determine
        # the the 'new name' for this FASTA file
        if [[ $bin_name == bin-unbinned ]]; then
            new_name="${new_culture_name}_UNBINNED.fa"
            echo "Working on $new_culture_name unbinned contigs (now called $new_name) ..."
        elif (( $(echo "$completion >= 50" | bc -l) && $(echo "$contamination < 5" | bc -l) )); then
            # did you realize that in the statement above we are using (( .. )) notation for our
            # conditional rather than [[ .. ]]? I know it is confusing, but the (( .. )) notation
            # is the shell's arithmetic construct and it is required to make numerical comparisons.
            # here is a useful summary of this confusing syntax:
            #
            #   https://superuser.com/questions/1533900/difference-between-and-or-and-in-bash
            #
            new_name=$(printf "${new_culture_name}_MAG_%03d.fa" $MAG_counter)
            echo "Working on $new_culture_name MAG $bin_name (now called $new_name) (C/R: $completion/$contamination) ..."
            MAG_counter=$(("$MAG_counter + 1"))
        else
            new_name=$(printf "${new_culture_name}_BIN_%03d.fa" $BIN_counter)
            echo "Working on $new_culture_name BIN $bin_name (now called $new_name) (C/R: $completion/$contamination) ..."
            BIN_counter=$(($BIN_counter + 1))
        fi

        # now the most tricky part -- we have a new name for our FASTA files,
        # but the CHECKM and GTDB-TAX files still have the old names. we need
        # to fix that, and sed is a great way to do it.
        sed -i '' "s/ms.*${bin_name}/$(basename $new_name .fa)/g" COMBINED-DATA/$new_culture_name-CHECKM.txt
        sed -i '' "s/ms.*${bin_name}/$(basename $new_name .fa)/g" COMBINED-DATA/$new_culture_name-GTDB-TAX.txt

        # finally we use anvi-script-reformat-fasta to both copy the files to the
        # right location, and rename individual sequences in them so each sequence
        # has a unique name -- having unique names for everything is a great idea
        # and will always have positive implications for downstream processes, but
        # if we didn't care about it, we could do this instead to copy files:
        #
        # cp $fasta_file COMBINED-DATA/$new_name
        #
        anvi-script-reformat-fasta $fasta_file \
                                   --simplify-names \
                                   --prefix $new_culture_name \
                                   --output-file COMBINED-DATA/$new_name \
                                   --quiet
    done
done
```
{% include CODEBLOCKFILENAME filename="generate-combined-data.sh" %}
