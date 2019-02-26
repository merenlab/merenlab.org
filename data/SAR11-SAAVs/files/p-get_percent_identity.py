'''loop over a bam file and get the edit distance to the reference genome stored in the NM tag
scale by aligned read length. works for bowtie2, maybe others. Adopted from
https://gigabaseorgigabyte.wordpress.com/2017/04/14/getting-the-edit-distance-from-a-bam-alignment-a-journey/'''

import sys
import pysam
import argparse
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

ap = argparse.ArgumentParser()

ap.add_argument("-b", required=False, help="bam filepaths comma separated (no spaces)")
ap.add_argument("-B", required=False, help="alternatively, a file of bam filepaths")
ap.add_argument("-p", required=False, action='store_true', help="provide if you want proper percent identity (unaligned bps are included in normalization)")
ap.add_argument("-u", required=False, action='store_true', help="if provided, histograms will be unnormalized")
ap.add_argument("-o", required=True, help="output filepath")
ap.add_argument("-r", required=False, default=None, help="If provided, only reads for specified references are considered. Comma separated.")
ap.add_argument("-R", required=False, default=None, help="If provided, only reads for specified references are considered. Filepath of list of references.")
ap.add_argument("-g", required=False, help="gene caller ids to report, comma-separated. requires -a flag")
ap.add_argument("-G", required=False, help="filepath of gene caller ids to report, single column file. requires -a flag")
ap.add_argument("-x", required=False, action='store_true', help="collapse histograms of individual genes into a single histogram per bam file. Valid only with -g and -G")
ap.add_argument("-a", required=False, help="output file of anvi-export-gene-calls, required for -g and -G, incompatible with -R")
ap.add_argument("-m", required=False, action='store_true', help="If provided, histogram will only be generated for the MEDIAN read length in each bam file. May want to use with -nummismatches")
ap.add_argument("-nummismatches", required=False, action='store_true', help="If provided, values are number of mismatches instead of percent identity. Highly recommended to use this with -m flag")
ap.add_argument("-mode", required=False, default='histogram', help="by default, this program reports histogram curves (-mode histogram). If ``-mode raw``, histograms are not computed and instead each read considered is a written with its percent identity value")

# These are parameters for mode = 'histogram'
ap.add_argument("-binsize", required=False, type=float, default=None, help="Size of each bin. Overrides behavior of -numbins")
ap.add_argument("-autobin", required=False, action='store_true', help="Size of each bin determined on a per bam basis (creates one bin for each mismatch, but you still must supply a -range value). -m is required for this mode")
ap.add_argument("-numbins", required=False, default=None, help="How many bins? default is 30")
ap.add_argument("-interpolate", default=None, required=False, help="How many points should form the interpolation and from where to where? Format is 'start,end,number'. (e.g. 67,100,200) If not provided, no interpolation. Required for autobin to normalize x-axis between bam files")
ap.add_argument("-range", required=False, default=None, help="What's the range? provide lower and upper bound comma-separated. default is 50,100")
ap.add_argument("-melted", required=False, action='store_true', help="Use melted output format. Required if using -autobin since each bam could have different bins")

# These are parameters for mode = 'raw'
ap.add_argument("-subsample", required=False, type=int, default=None, help="how many reads do you want to subsample? You probably don't need more than 50,000.")

args = vars(ap.parse_args())
sc = lambda parameters: any([args.get(parameter, False) for parameter in parameters])

raw_mode_parameters = ['subsample']
histogram_mode_parameters = ['binsize', 'autobin', 'numbins', 'interpolate', 'range', 'melted']
if args.get("mode") == 'histogram':
    if sc(raw_mode_parameters):
        raise Exception(" you are using mode = histogram. Don't use these parameters: {}".format(raw_mode_parameters))

    # defaults
    if not args.get("numbins"):
        args['numbins'] = 30
    if not args.get("range"):
        args['range'] = "50,100"

if args.get("mode") == 'raw':
    if sc(histogram_mode_parameters):
        raise Exception(" you are using mode = raw. Don't use these parameters: {}".format(histogram_mode_parameters))



# checks
if args.get("g") and args.get("G"):
    raise Exception("use -g or -G")
if (args.get("g") or args.get("G")) and not args.get("a"):
    raise Exception("provide -a")
if (args.get("g") or args.get("G")) and (args.get("R") or args.get("r")):
    raise Exception("no point providing -R/-r if using gene calls")
if args.get("x") and not (args.get('g') or args.get('G')):
    raise Exception("you need to specify -g or -G to use -x")
if args.get("b") and args.get('B'):
    raise Exception("specify either b or B, not both")
if not args.get("b") and not args.get('B'):
    raise Exception("Specify one of b or B")
if args.get('autobin') and not args.get("melted"):
    raise Exception("You can't autobin without using -melted output format.")
if args.get('autobin') and not args.get("m"):
    raise Exception("You can't autobin without using -m.")
if args.get('autobin') and not args.get("interpolate"):
    raise Exception("You can't autobin without using -interpolate.")

if args.get("g"):
    genes = [int(x) for x in args['g'].split(',')]
elif args.get("G"):
    print(args['G'])
    genes = [int(x.strip()) for x in open(args['G']).readlines()]
else:
    genes = None

if args.get("a"):
    gene_info = pd.read_csv(args['a'], sep='\t')
    if genes:
        gene_info = gene_info[gene_info['gene_callers_id'].isin(genes)]
else:
    gene_info = None

if args.get('b'):
    bam_filepaths = args['b'].split(",") # comma separated bam file paths
    if bam_filepaths[-1] == '':
        bam_filepaths = bam_filepaths[:-1]
elif args.get('B'):
    bam_filepaths = [x.strip() for x in open(args['B']).readlines()]
    if not bam_filepaths:
        print('no bam files provided. nothing to do.')
        sys.exit()

proper_pident = args['p']
output = args['o']
if args.get('r'):
    reference_names = args['r'].split(",") # comma separated reference names
elif args.get('R'):
    reference_names = [x.strip() for x in open(args['R']).readlines()]
else:
    reference_names = [None]
normalize = True if not args['u'] else False

############################################################################################################################

# preamble
attr = 'query_alignment_length' if not proper_pident else 'query_length'

if args.get('mode') == 'histogram':
    range_low, range_hi = [int(x) for x in args['range'].split(",")]

    if not args.get("autobin"):
        if not args.get("binsize"):
            bins_template = np.linspace(range_low, range_hi, int(args['numbins']), endpoint=True)
        else:
            bins_template = np.arange(range_hi, range_low - float(args['binsize']), -float(args['binsize']))[::-1]
        if not args.get("nummismatches"):
            bins_shifted = bins_template[1:] # include 100% as a datapoint
        else:
            bins_shifted = bins_template[:-1] # include 0 mismatches as a datapoint

    if args.get("interpolate"):
        interp_low, interp_hi, interp_num = [int(x) for x in args['interpolate'].split(",")]
        interp_low += 0.5
        bins_interp = np.linspace(interp_low, interp_hi, interp_num)
    else:
        bins_interp = np.array([])

I = lambda bins_shifted, counts, bins_interp: interp1d(bins_shifted, counts, kind='cubic')(bins_interp) if args['interpolate'] else counts
J = lambda true_or_false, length: length == median_length if true_or_false else True
K = lambda true_or_false, read: read.get_tag("NM") if true_or_false else 100*(1 - float(read.get_tag("NM")) / read.__getattribute__(attr))

if args.get("mode") == 'histogram':
    percent_identity_hist = {'value':[], 'percent_identity':[], 'id':[]} if args.get("melted") else {}
if args.get("mode") == 'raw':
    percent_identity_hist = {'value':[], 'id':[]}

for bam in bam_filepaths:
    print('working on {}...'.format(bam))
    bam_name = bam.split("/")[-1].replace(".bam", "")
    samfile = pysam.AlignmentFile(bam, "rb")

    if args.get("m"):
        i = 0
        read_lengths = []
        for read in samfile.fetch():
            read_lengths.append(read.query_length)
            i += 1
        array = np.array(read_lengths)
        median_length = int(np.median(array))
        second_median = int(np.median(array[array != float(median_length)]))
        third_median = int(np.median(array[(array != float(median_length)) & (array != float(second_median))]))
        print('median length was {}, second was {}, third was {}'.format(median_length, second_median, third_median))

        # only if autobinning
        if args.get("autobin"):
            if not args.get("nummismatches"):
                binsize = 100*(1. / median_length)
                bins_template = np.arange(range_hi, range_low - binsize, -binsize)[::-1]
                bins_template += binsize * 1e-4
                bins_template[-1] = 100
                bins_shifted = bins_template[1:] # include 100% as a datapoint
            else:
                bins_template = np.arange(range_hi, range_low - 1, -1)[::-1]
                bins_shifted = bins_template[:-1] # include 0 mismatches as a datapoint

    if gene_info is not None:
        if not args.get('x'):
            # each gene gets its own histogram
            for index, row in gene_info.iterrows():
                percent_identities = np.array([K(args.get("nummismatches"), read) for read in samfile.fetch(row['contig'], int(row['start']), int(row['stop'])) if J(args.get("m"), read.query_length)])
                id_name = bam_name + "_" + str(row['gene_callers_id'])
                if args.get("mode") == 'histogram':
                    counts, _ = np.histogram(percent_identities, bins=bins_template, density=normalize)
                    if args.get("melted"):
                        value = list(I(bins_shifted, counts, bins_interp))
                        percent_identity = list(bins_shifted if not args.get("interpolate") else bins_interp)
                        percent_identity_hist['id'].extend([id_name] * len(value))
                        percent_identity_hist['value'].extend(value)
                        percent_identity_hist['percent_identity'].extend(percent_identity)
                    else:
                        percent_identity_hist[id_name] = I(bins_shifted, counts, bins_interp)
                elif args.get("mode") == 'raw':
                    if args.get("subsample") and args['subsample'] < len(percent_identities):
                        percent_identities = np.random.choice(percent_identities, args.get('subsample'), replace=False)
                    value = percent_identities.tolist()
                    percent_identity_hist['id'].extend([id_name] * len(value))
                    percent_identity_hist['value'].extend(value)
        else:
            # histograms for each gene are collapsed into a single histogram
            percent_identities = []
            for index, row in gene_info.iterrows():
                percent_identities.extend([K(args.get("nummismatches"), read) for read in samfile.fetch(row['contig'], int(row['start']), int(row['stop'])) if J(args.get("m"), read.query_length)])
            percent_identities = np.array(percent_identities)
            id_name = bam_name
            if args.get("mode") == 'histogram':
                counts, _ = np.histogram(percent_identities, bins=bins_template, density=normalize)
                if args.get("melted"):
                    value = list(I(bins_shifted, counts, bins_interp))
                    percent_identity = list(bins_shifted if not args.get("interpolate") else bins_interp)
                    percent_identity_hist['id'].extend([id_name] * len(value))
                    percent_identity_hist['value'].extend(value)
                    percent_identity_hist['percent_identity'].extend(percent_identity)
                else:
                    percent_identity_hist[id_name] = I(bins_shifted, counts, bins_interp)
            elif args.get("mode") == 'raw':
                if args.get("subsample") and args['subsample'] < len(percent_identities):
                    percent_identities = np.random.choice(percent_identities, args.get('subsample'), replace=False)
                value = percent_identities.tolist()
                percent_identity_hist['id'].extend([id_name] * len(value))
                percent_identity_hist['value'].extend(value)
    else:
        percent_identities = []
        for reference_name in reference_names:
            # 1 minus the ratio of mismatches to the length of the read/alignment
            percent_identities.extend([K(args.get("nummismatches"), read) for read in samfile.fetch(reference=reference_name) if J(args.get("m"), read.query_length)])
        percent_identities = np.array(percent_identities)
        print("{} reads considered".format(len(percent_identities)))

        # create a histogram
        id_name = bam_name
        if args.get("mode") == 'histogram':
            counts, _ = np.histogram(percent_identities, bins=bins_template, density=normalize)
            if args.get("melted"):
                value = list(I(bins_shifted, counts, bins_interp))
                percent_identity = list(bins_shifted if not args.get("interpolate") else bins_interp)
                percent_identity_hist['id'].extend([id_name] * len(value))
                percent_identity_hist['value'].extend(value)
                percent_identity_hist['percent_identity'].extend(percent_identity)
            else:
                percent_identity_hist[id_name] = I(bins_shifted, counts, bins_interp)
        if args.get("mode") == 'raw':
                if args.get("subsample") and args['subsample'] < len(percent_identities):
                    percent_identities = np.random.choice(percent_identities, args.get('subsample'), replace=False)
                value = percent_identities.tolist()
                percent_identity_hist['id'].extend([id_name] * len(value))
                percent_identity_hist['value'].extend(value)
    samfile.close()
    print("")

# save the file
percent_identity_hist = pd.DataFrame(percent_identity_hist).reset_index(drop=True)
if not args.get("melted") and args.get("mode") == 'histogram':
    percent_identity_hist['percent_identity' if not args['nummismatches'] else 'number_of_mismatches'] = bins_interp if args['interpolate'] else bins_shifted
percent_identity_hist.to_csv(output, sep="\t", index=False)

