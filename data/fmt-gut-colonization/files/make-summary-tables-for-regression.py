import pandas as pd

# Get group non-specific dfs:
transplants_raw = pd.read_csv('metadata-transplants.txt', sep='\t', index_col=False)

for group in ['DA', 'DB']:

    # Get group specific dfs:
    colo_df = pd.read_csv(f'colonized-{group}.txt', sep='\t', index_col=None)
    no_colo_df = pd.read_csv(f'did-not-colonize-{group}.txt', sep='\t', index_col=None)
    transplants = transplants_raw.loc[transplants_raw['Group'] == group]
    countries_summary = pd.read_csv(f'detection-global-by-country-{group}.txt', sep='\t', index_col=0)
    coverage_raw = pd.read_csv(f'mean-cov-{group}.txt', sep='\t', index_col=None)
    coverage = coverage_raw.melt(id_vars=['bins'])

    # Combine colo and no_colo data together:
    colo_df['outcome'] = 'colonization'
    no_colo_df['outcome'] = 'no_colonization'
    outcomes = pd.concat([no_colo_df, colo_df])

    # Merge dfs together:
    merge1 = pd.merge(outcomes, transplants[['Recipient', 'Sample Name', 'FMT Method']], left_on='recipient', right_on='Recipient', how='left')
    merge1.drop('Recipient', axis=1, inplace=True)
    merge1.rename(columns={'Sample Name':'transplant_sample', 'FMT Method':'fmt_method'}, inplace=True)
    merge2 = pd.merge(merge1, coverage, left_on=['MAG', 'transplant_sample'], right_on=['bins', 'variable'], how='left')
    merge2.rename(columns={'value':'transplant_mean_cov_Q2Q3'}, inplace=True)
    merge2.drop(columns=['variable', 'bins'], inplace=True)
    merge3 = pd.merge(merge2, countries_summary, left_on='MAG', right_on='bins', how='left')

    # Save:
    merge3.to_csv(f'summary-for-regression-{group}.txt', sep='\t', index=False)
