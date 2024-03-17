#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse
import numpy as np

from itertools import chain
from collections import Counter

import anvio.utils as utils
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

regions_of_interests = {'repA': (1809, 2467), 'mobA': (300, 1404)}
min_avearge_coverage = 1

def main(args):
    input_file_paths = args.input_files

    regions_max = max([max(e) for e in regions_of_interests.values()])
    regions_min = min([min(e) for e in regions_of_interests.values()])

    d = {}
    for input_file_path in input_file_paths:
        filesnpaths.is_file_exists(input_file_path)
        filesnpaths.is_file_tab_delimited(input_file_path, expected_number_of_fields=5)

        # alternative input to consider. will appear in `d` later:
        alt = input_file_path.split('.txt')[0]

        contents_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_path)

        for e in contents_dict.values():
            pos = int(e['nt_position'])
            if pos > regions_min and pos < regions_max:
                # this is within the min/max range of thins we are interested but
                # we are not sure if it is relevant yet
                for region in regions_of_interests:
                    if pos > min(regions_of_interests[region]) and pos < max(regions_of_interests[region]):
                        # we are in 'region'
                        sample = e['sample_name']

                        if sample not in d:
                            d[sample] = {}

                        if alt not in d[sample]:
                            d[sample][alt] = {}

                        if region not in d[sample][alt]:
                            d[sample][alt][region] = []

                        d[sample][alt][region].append(int(e['coverage']))
                    else:
                        continue

    alternatives = {}
    mean_coverages = {}
    for sample in d:
        mean_coverages[sample] = Counter()
        for alt in d[sample]:

            # store the average coverage of all regions for the alternative
            mean_coverages[sample][alt] = np.mean(list(chain(*d[sample][alt].values())))

            for region in d[sample][alt]:
                d[sample][alt][region] = {'cov': np.mean(d[sample][alt][region]), 'det': len([p for p in d[sample][alt][region] if p > 1]) / len(d[sample][alt][region])}


        # remove alternatives if they don't have prper mean coverage
        alternatives_to_remove_due_to_min_coverage = [alt for alt in mean_coverages[sample] if mean_coverages[sample][alt] < min_avearge_coverage]
        for alt in alternatives_to_remove_due_to_min_coverage:
            d[sample].pop(alt)

        if not len(d[sample]):
            alternatives[sample] = {'plasmid_version': 'NONE_DUE_TO_MIN_COVERAGE'}
            continue

        regions_for_consideration = []
        # first check for detection. if any of regions of interest has detection less than 1.0,
        # they are out of consideration:
        for alt in d[sample]:
            if len([r for r in d[sample][alt] if d[sample][alt][r]['det'] > 0.9]) == len(regions_of_interests):
                regions_for_consideration.append(region)

        if not len(regions_for_consideration):
            alternatives[sample] = {'plasmid_version': 'NONE_DUE_TO_NO_REG_DETECTION'}
            continue

        # for the remaining ones, we want to compare the ratio between alternatives
        for alt in d[sample]:
            d[sample][alt]['ratios'] = d[sample][alt]['mobA']['cov'] / d[sample][alt]['repA']['cov']

        ratios = [(alt, d[sample][alt]['ratios']) for alt in d[sample]]
        ratios_within_range = [(r[0], abs(1 - r[1])) for r in ratios if r[1] > 0.1 and r[1] < 4]

        if not len(ratios_within_range):
            alternatives[sample] = {'plasmid_version': 'NONE_DUE_TO_BAD_RATIOS'}
            continue

        # alternative version of the plasmid for this metagenome is determined here:
        alternative = sorted(ratios_within_range, key = lambda x: x[1])[0][0]
        alternatives[sample] = {'plasmid_version': alternative}

    utils.store_dict_as_TAB_delimited_file(alternatives, 'ALTERNATIVES.txt', headers=['metagenome', 'plasmid_version'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process multiple split nt coverages files')

    parser.add_argument('input_files', metavar = 'SPLIT_NT_COVERAGE_FILES', nargs='+',
                        help = "Split nt coverages to process")

    args = parser.parse_args()

    try:
        main(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)
