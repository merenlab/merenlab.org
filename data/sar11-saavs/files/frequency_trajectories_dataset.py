import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import argparse
import sys
import os

ap = argparse.ArgumentParser()
ap.add_argument('--metadata', '-m', default='TARA_metadata.txt')
ap.add_argument('--variability', '-vv', required=True)
ap.add_argument('--output', '-o', required=True)
ap.add_argument('--genes', '-g', default=None)
ap.add_argument('--positions', '-p', default=None)
ap.add_argument('--full', '-f', action='store_true')
ap.add_argument('--plot', action='store_true')
ap.add_argument('--correlation-variable', default='Temperature')
args = ap.parse_args()

if args.plot and not args.full:
    print('must use --full if you want to plot')
    sys.exit()
if args.plot and not (args.genes or args.positions):
    print('you cant plot everything. Pick genes or positions')
    sys.exit()
if args.genes and args.positions:
    print('pick genes or positions, not both')

meta = pd.read_csv(args.metadata, sep='\t').rename(columns={'Temperature (deg C)': 'Temperature', 'Sample Id': 'sample_id'})
df = pd.read_csv(args.variability, sep='\t')
df['Temperature'] = df['sample_id'].map(dict(zip(meta['sample_id'], meta['Temperature'])))
df['NO2NO3'] = df['sample_id'].map(dict(zip(meta['sample_id'], meta['NO2NO3 (umol/L)'])))
df['Nitrates'] = df['sample_id'].map(dict(zip(meta['sample_id'], meta['Nitrates (umol/L)'])))
df['Depth'] = df['sample_id'].map(dict(zip(meta['sample_id'], meta['Depth Id'])))
number_of_samples = df['sample_id'].nunique()

if args.genes:
    genes = [x.strip() for x in open(args.genes).readlines()]
    df = df[df['corresponding_gene_call'].isin(genes)]
    positions = list(df['unique_pos_identifier'].unique())
elif args.positions:
    positions = [x.strip() for x in open(args.positions).readlines()]
else:
    print("you didn't provide it genes or a positions file. Don't try and visualize this")
    positions = list(df['unique_pos_identifier'].unique())
genes = list(df.loc[df['unique_pos_identifier'].isin(positions), 'corresponding_gene_call'].unique())
df = df[df['corresponding_gene_call'].isin(genes)]

convert = [('Ala', 'A'), ('Arg', 'R'), ('Asn', 'N'), ('Asp', 'D'), ('Cys', 'C'), ('Gln', 'Q'), ('Glu', 'E'), ('Gly', 'G'), ('His', 'H'), ('Ile', 'I'), ('Leu', 'L'), ('Lys', 'K'), ('Met', 'M'), ('Phe', 'F'), ('Pro', 'P'), ('Ser', 'S'), ('Thr', 'T'), ('Trp', 'W'), ('Tyr', 'Y'), ('Val', 'V')]
aas = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
groups = {"Positive": ["R", "H", "K"],
          "Negative": ["D", "E"],
          "Polar": ["S", "T", "N", "Q"],
          "Special": ["C", "G", "P"],
          "Hydrophobic": ["A", "V", "M", "F", "I", "L", "W", "Y"]}
groups = {{val: key for key, val in convert}[aa]: group \
                    for group, aas in groups.items() for aa in aas}

if args.plot:
    try: 
        os.mkdir('aast_movement')
    except:
        pass

output = {'rvalue'                : [],
          'R2'                    : [],
          'pvalue'                : [],
          'slope'                 : [],
          'absolute_slope'        : [],
          'aa1'                   : [],
          'aa2'                   : [],
          'AAST'                  : [],
          'class'                 : [],
          'codon'                 : [],
          'gene'                  : [],
          'unique_pos_identifier' : []}
if args.full:
    output.update({
        'frequency1'  : [],
        'frequency2'  : [],
        'coverage'    : [],
        'Temperature' : [],
        'NO2NO3' : [],
        'Nitrates' : [],
        'Depth' : [],
    })
if 'solvent_acc' in df.columns:
    output['solvent_accessibility'] = []

number_of_aast_positions = len(positions)
print('total number of positions considered: {}'.format(number_of_aast_positions))
ind = 0
for gene, df_subset in df.groupby('corresponding_gene_call'):
    gene_positions = df_subset.loc[df_subset['unique_pos_identifier'].isin(positions), 'unique_pos_identifier'].unique()
    for AAST_position in gene_positions:
        temp = df_subset[df_subset['unique_pos_identifier'] == AAST_position].sort_values(by='Temperature')
        codon = temp['codon_order_in_gene'].iloc[0]
        if 'solvent_acc' in df.columns:
            solvent_accessibility = temp['solvent_acc'].iloc[0]

        #  try to pick the most common. If the most common is an amino acid verses itself, pick the
        #  second most common, and so on and so forth. Also demand that the stop codon is not one of
        #  them. If you can't find an AAST meeting this criteria, move onto the next position
        AASTs = temp['competing_aas'].value_counts()
        found = False
        while True:
            if AASTs.index[0][:3] != AASTs.index[0][3:] and 'STP' not in AASTs.index[0]:
                found = True
                AAST = AASTs.index[0]
                break
            else:
                if len(AASTs) == 1:
                    break
                AASTs = AASTs.iloc[1:]
        if not found:
            continue

        aa1 = AAST[:3]
        aa2 = AAST[3:]

        position = temp['unique_pos_identifier'].iloc[0]

        temperature = temp['Temperature']
        depth = temp['Depth']
        NO2NO3 = temp['NO2NO3']
        nitrates = temp['Nitrates']
        aa1_total = temp[aa1] + temp[aa2]
        aa1_freq = temp[aa1] / aa1_total
        aa2_freq = temp[aa2] / aa1_total

        # if a metagenome does not contain coverage for both of the amino acids, the whole position
        # is discarded
        no_competition = np.isnan(aa1_freq)
        aa1_freq = aa1_freq[~no_competition]
        aa2_freq = aa2_freq[~no_competition]
        aa1_total = aa1_total[~no_competition]
        temperature = temperature[~no_competition]
        NO2NO3 = temperature[~no_competition]
        nitrates = temperature[~no_competition]
        depth = depth[~no_competition]
        if len(aa1_freq) < number_of_samples:
            continue

        if args.correlation_variable == 'NO2NO3':
            linear_regression = linregress(NO2NO3, aa1_freq)
        elif args.correlation_variable == 'Nitrates':
            linear_regression = linregress(nitrates, aa1_freq)
        elif args.correlation_variable == 'Temperature':
            linear_regression = linregress(temperature, aa1_freq)
        else:
            print('not a valid correlation variable')
            sys.exit()

        if args.plot:
            plt.scatter(temperature, aa1_freq, label=aa1)
            plt.scatter(temperature, aa2_freq, label=aa2)
            plt.legend(loc='best')
            plt.savefig('aast_movement/{}_{}.png'.format(gene, codon))
            plt.close()

        if args.full:
            output['Temperature'].extend(temperature)
            output['Depth'].extend(depth)
            output['NO2NO3'].extend(NO2NO3)
            output['Nitrates'].extend(nitrates)
            output['frequency1'].extend(aa1_freq)
            output['frequency2'].extend(aa2_freq)
            output['coverage'].extend(aa1_total)
            output['rvalue'].extend([linear_regression.rvalue] * q)
            output['R2'].extend([linear_regression.rvalue**2] * q)
            output['pvalue'].extend([linear_regression.pvalue] * q)
            output['slope'].extend([linear_regression.slope] * q)
            output['absolute_slope'].extend([abs(linear_regression.slope)] * q)
            output['aa1'].extend([aa1] * q)
            output['aa2'].extend([aa2] * q)
            output['AAST'].extend([aa1+aa2] * q)
            output['class'].extend([groups[aa1]+groups[aa2]] * q)
            output['codon'].extend([codon] * q)
            output['gene'].extend([gene] * q)
            output['unique_pos_identifier'].extend([AAST_position] * q)
            if 'solvent_acc' in df.columns:
                output['solvent_accessibility'].extend([solvent_accessibility] * q)
        else:
            output['rvalue'].append(linear_regression.rvalue)
            output['R2'].append(linear_regression.rvalue**2)
            output['pvalue'].append(linear_regression.pvalue)
            output['slope'].append(linear_regression.slope)
            output['absolute_slope'].append(abs(linear_regression.slope))
            output['aa1'].append(aa1)
            output['aa2'].append(aa2)
            output['AAST'].append(aa1+aa2)
            output['class'].append(groups[aa1]+groups[aa2])
            output['codon'].append(codon)
            output['gene'].append(gene)
            output['unique_pos_identifier'].append(AAST_position)
            if 'solvent_acc' in df.columns:
                output['solvent_accessibility'].append(solvent_accessibility)

        ind += 1

output = pd.DataFrame(output)
output.to_csv(args.output, sep='\t', index=False)
