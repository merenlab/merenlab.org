---
layout: page
noleftpanel: true
title: "Solution for EXC-003"
author: meren
date: February, 2025
geometry: "left=3cm,right=3cm,top=2cm,bottom=2cm"
colorlinks: true
urlcolor: blue
monofont: DejaVuSansMono.ttf
---

[Go back to the question](../../#exc-003).

```sh
## A simple shell script to learn about the properties of sequences
## in a given FASTA file

set -e

## Check if the user provided a FASTA file
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <FASTA file>"
    exit 1
fi

FASTA_FILE="$1"

## Check if the file exists and is readable
if [[ ! -f "$FASTA_FILE" || ! -r "$FASTA_FILE" ]]; then
    echo "Error: File does not exist or is not readable."
    exit 1
fi

## Validate if it's a FASTA file (basic check for '>' at the start)
if ! grep -q "^>" "$FASTA_FILE"; then
    echo "Error: File does not appear to be a FASTA file."
    exit 1
fi

## Process the FASTA file with AWK
awk '
BEGIN {
    seq_count = 0;
    total_length = 0;
    max_length = 0;
    min_length = 99999999;
    gc_count = 0;
    sequence = "";
}

## Identify sequence headers
/^>/ {
    if (sequence != "") {
        seq_length = length(sequence);
        total_length += seq_length;
        if (seq_length > max_length) max_length = seq_length;
        if (seq_length < min_length) min_length = seq_length;

        # Count GC content
        gc_count += gsub(/[GgCc]/, "", sequence);
    }
    seq_count++;
    sequence = "";
    next;
}

## Append sequence lines
{
    sequence = sequence $0;
}

END {
    if (sequence != "") {
        seq_length = length(sequence);
        total_length += seq_length;
        if (seq_length > max_length) max_length = seq_length;
        if (seq_length < min_length) min_length = seq_length;

        # Count GC content
        gc_count += gsub(/[GgCc]/, "", sequence);
    }

    # Calculate average sequence length
    avg_length = (seq_count > 0) ? total_length / seq_count : 0;

    # Compute GC percentage
    gc_percent = (total_length > 0) ? (gc_count / total_length) * 100 : 0;

    # Print statistics
    print "FASTA File Statistics:";
    print "----------------------";
    print "Number of sequences: " seq_count;
    print "Total length of sequences: " total_length;
    print "Length of the longest sequence: " max_length;
    print "Length of the shortest sequence: " min_length;
    print "Average sequence length: " avg_length;
    print "GC Content (%): " gc_percent;
}' "$FASTA_FILE"

```
{% include CODEBLOCKFILENAME filename="fasta-file-processor.sh" %}
