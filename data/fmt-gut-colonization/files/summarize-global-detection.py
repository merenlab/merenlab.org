#!/usr/bin/env python3
import pandas as pd

# A function to summarize detection of bins in different groups designated by 3-letter sample prefixes

def summarize(detection_file, output):
    # Read detection table, with bin names as indexes:
    detection = pd.read_csv(detection_file, sep='\t', index_col='bins')
    
    # Set detection values less than 0.25 to 0 and
    # set detection values greater than or equal to 0.25 to 1:
    detection[detection < 0.25] = 0
    detection[detection >= 0.25] = 1
    
    # Make list of all samples, alphabetically just because:
    samples = sorted(detection.columns.values.tolist())
    
    # Make a set of the first 3 letters of all samples:
    prefixes = [sample[:3] for sample in samples]
        
    # Turn the list of prefixes into a set so it only contains unique values:
    groups = sorted(set(prefixes))

    # Make a dictionary of how many samples are in each group:
    group_counts = {}
    for group in groups:
        group_counts.update({group: prefixes.count(group)})
    
    # Make a summary df with the same indexes as the detection df:
    summary = pd.DataFrame(index=detection.index)

    # Add new columns stating 1) how many times each bin is detected within a group
    # and 2) the fraction of samples in a group containing that bin:
    for group in groups:
        group_columns = detection.filter(regex=group).copy()

        summary[f'{group}_count_detec'] = group_columns.sum(axis=1)

        summary[f'{group}_frac_detec'] = group_columns.sum(axis=1)/len(group_columns.columns)
        
    # Save final summary table as tsv:
    summary.to_csv(output, sep='\t', index=True) 

## Using this function to summarize the detection of the Canada FMT CDI IBD donor-derived MAGs in healthy adult gut metagenomes from different countries

summarize('detection-global-DA.txt', 'detection-global-by-country-DA.txt')

summarize('detection-global-DB.txt', 'detection-global-by-country-DB.txt')
