#!/usr/bin/env python3
import pandas as pd

for group in ['A', 'B']:
    mags_summary = pd.read_csv(f"FMT_DONOR_{group}_AND_RECIPIENTS/SUMMARY/bins_summary.txt", sep="\t", index_col=0)
    mags = list(mags_summary.index)
    hmm_means = {}

    for mag in mags:
        # Get gene caller ids for Campell et al. HMMs:
        mag_path = f"FMT_DONOR_{group}_AND_RECIPIENTS/SUMMARY/bin_by_bin/{mag}/{mag}-"
        mag_hmm_seqs = f"{mag_path}Campbell_et_al-hmm-sequences.txt"
        mag_hmm_ids = []

        with open(mag_hmm_seqs, "r") as hmms_file:
            for line in hmms_file:
                stripped_line = line.strip()
                gene_caller_id_start = stripped_line.find("gene_callers_id:") + 16

                if gene_caller_id_start != 15:
                    gene_caller_id_stop = stripped_line.find("|", gene_caller_id_start)
                    gene_caller_id = stripped_line[gene_caller_id_start:gene_caller_id_stop]
                    mag_hmm_ids.append(int(gene_caller_id))

        # Get coverages of HMM gene caller ids in each sample:
        mag_gene_covs = pd.read_csv(f"{mag_path}gene_non_outlier_coverages.txt", sep='\t', index_col=0)

        if not hmm_means:
            samples = mag_gene_covs.columns.tolist()

        mag_hmm_covs = mag_gene_covs.loc[ mag_hmm_ids , : ]
        mag_hmm_covs.loc['mean'] = mag_hmm_covs.mean()
        mag_hmm_means = mag_hmm_covs.loc['mean'].tolist()
        hmm_means[mag] = mag_hmm_means

    hmm_means_df = pd.DataFrame(hmm_means,index=samples).T
    hmm_means_df.to_csv(f"scg-cov-D{group}.txt", sep='\t', index=True)
