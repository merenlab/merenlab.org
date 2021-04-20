import pandas as pd

taxon_levels = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')

groups_to_merge = ('DA', 'DB', 'CAN')

group_name = "_".join(groups_to_merge)

for level in taxon_levels:

    dfs = []

    for group in groups_to_merge:
        df = pd.read_csv(f'taxonomy-tables/{group}_t_{level}.txt', sep='\t', index_col='key')
        dfs.append(df)

    final = pd.concat(dfs, axis=0, ignore_index=False, sort=True)
    final.to_csv(f'{group_name}_t_{level}.txt', sep='\t', na_rep=0)
