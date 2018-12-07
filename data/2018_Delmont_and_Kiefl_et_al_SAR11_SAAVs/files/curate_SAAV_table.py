import pandas as pd

saav_table_fname = "S-LLPA_SAAVs_20x_10percent_departure.txt"
nonoutlier_gene_cov_fname = "ALL_SPLITS-gene_non_outlier_coverages.txt"

# load table
saav_table = pd.read_csv(saav_table_fname, sep='\t', header=0, index_col=False)

# neglect inter-gene stop codons, which have invalid BLOSUM scores
saav_table.dropna(inplace=True, axis=0)

# filter all entries whose departure from consensus <= 0.1
saav_table = saav_table[saav_table["departure_from_consensus"] > 0.1]

# append the average non-outlier gene coverage information to each gene as a column in SAAV table
avg_gene_cov = pd.read_csv(nonoutlier_gene_cov_fname, sep="\t", index_col=0)
avg_gene_cov["corresponding_gene_call"] = avg_gene_cov.index
avg_gene_cov = pd.melt(avg_gene_cov, id_vars="corresponding_gene_call", var_name="sample_id", value_name="nonoutlier_gene_cov")
saav_table = saav_table.merge(avg_gene_cov)

# add rel_diff_from_mean_gene_cov
saav_table["rel_diff_from_mean_gene_cov"] = (saav_table.coverage - saav_table.nonoutlier_gene_cov) / saav_table.nonoutlier_gene_cov

# simplify sample names
saav_table["sample_id"] = saav_table["sample_id"].str.replace("_BOWTIE2","").str.replace("SAR11_","")

# save table
saav_table.to_csv("S-LLPA_SAAVs_20x_10percent_departure_curated.txt", sep='\t', index=False)
