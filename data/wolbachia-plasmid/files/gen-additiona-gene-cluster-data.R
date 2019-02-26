#!/usr/bin/env Rscript

library("dplyr")

# read the gene clusters summary table which is the output of
# anvi-summarize command run prior:
gene_clusters_df <- read.table(file = 'Wolbachia-Pan/SUMMARY/Wolbachia-Pan_gene_clusters_summary.txt',
                               header = TRUE,
                               sep = "\t",
                               quote = "")

# remove singleton gene clusters unique to wPip since they have no coverage values
gene_clusters_df_of_wPip <- gene_clusters_df[gene_clusters_df$genome_name == "wPip", ]

# read the blast output for best WO hits among wPip genes
WO_in_wPip_df <- read.table(file = 'WO_in_wPip_best_hit.txt',
                            header = TRUE, sep = "\t", quote = "")

# learn which gene cluster ids corrpespond to wPip gene calls
gene_clusters_df_of_wPip_phages=left_join(gene_clusters_df_of_wPip,
                                          WO_in_wPip_df,
                                          by="gene_callers_id")

# keep only columns of interest
gene_clusters_df_of_wPip_phages_short <- gene_clusters_df_of_wPip_phages[, colnames(gene_clusters_df_of_wPip_phages) %in% c('WO_assignment', 'pct_alignment', 'gene_cluster_id')]

# now we need to generate coverage values for each gene cluster.
# O07 had the most number of gene clusters, presenting a meaningful
# target to have coverage values for as many gene clusters as possible.
# here, we first get the gene clusters for O07
gene_clusters_df_of_O07 <- gene_clusters_df[gene_clusters_df$genome_name == "Wolbachia_O07", ]

# then remove all columns but two
gene_clusters_df_of_O07 <- gene_clusters_df_of_O07[, colnames(gene_clusters_df_of_O07) %in% c('gene_cluster_id', 'gene_callers_id')]

# learn gene coverages for all four Wolbachia genomes in the O07 metagenome using
# the summary of the merged profile database for O07:
gene_coverages_df_of_O07 <- read.table(file='06_MERGED/Culex_O07/SUMMARY/bin_by_bin/Wolbachia/Wolbachia-gene_coverages.txt',
                                       header=TRUE, sep="\t", quote="")

# anvi'o adds ribosomal RNA genes ad hoc after HMM runs, so they appear among gene calls but not in the
# pangenome which only uses ORFs identified through prodigal. we will remove them bofore continuing:
missing_genes <- setdiff(gene_coverages_df_of_O07$gene_callers_id,
                         gene_clusters_df_of_O07$gene_callers_id)
gene_coverages_df_of_O07 <- gene_coverages_df_of_O07[!(gene_coverages_df_of_O07$gene_callers_id %in% missing_genes), ]

# now we can merge the two;
gene_clusters_df_of_O07_with_covs=full_join(gene_clusters_df_of_O07,
                                            gene_coverages_df_of_O07,
                                            by="gene_callers_id")

# the following takes the mean coverage of genes per gene cluster
Mean_coverage_gene_clusters_df_of_O07 <- gene_clusters_df_of_O07_with_covs %>%
                                         group_by(gene_cluster_id) %>%
                                         summarize(mean(Wolbachia_O03),
                                                   mean(Wolbachia_O07),
                                                   mean(Wolbachia_O11),
                                                   mean(Wolbachia_O12))

# Join WO assignment and Coverage so that we get GC common to wPip and
# O07 with 4 additional layers
gene_clusters_df <- full_join(gene_clusters_df_of_wPip_phages_short,
                              Mean_coverage_gene_clusters_df_of_O07,
                              by="gene_cluster_id")

# give better names for each column
names(gene_clusters_df) <- c('gene_cluster_id',
                             'WO_assignment',
                             'pct_alignment',
                             'cov_in_003',
                             'cov_in_007',
                             'cov_in_011',
                             'cov_in_012')

# Here we have a data frame that lists for each gene cluster (1) the mean coverage values
# per gene per genome, (2) WO assignments if genes match to wPip prophages, and (3) the percent
# alignment of the query gene seqeunce so we can filter out weak hits. the following
# section of the code works on this data frame to produce a final additional item
# data file to be incorporated into the pan genome.

# a better max function that the default max function. just ignore:
better_max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

output_df <- data.frame(gene_cluster_id=character(),
                        WO_assignmnet=character(),
                        pct_alignment=numeric(),
                        cov_in_003=numeric(),
                        cov_in_007=numeric(),
                        cov_in_011=numeric(),
                        cov_in_012=numeric(),
                        stringsAsFactors=FALSE)

# here we go through every gene cluster, pick the most frequent WO assignment if there are any,
# and then discard the ones resulted from hits where less than 90% of the query sequence
# were aligned
N = 1
for (gene_cluster in levels(gene_clusters_df$gene_cluster_id)){
  # gene_cluster <- "GC_00001082"
  df <- gene_clusters_df[gene_clusters_df$gene_cluster_id == gene_cluster, ]
  freqs <- as.data.frame(table(df$WO_assignment))
  freqs <- freqs[freqs$Freq > 0, ]
  most_frequent_assignment <- as.character(freqs[order(c(-freqs$Freq)), ]$Var1[1])

  best_pct_alignment_score <- better_max(df$pct_alignment)

  # remove weak WO assignments
  if (!is.na(df$pct_alignment[1])){
    if (df$pct_alignment[1] < 90){
      most_frequent_assignment <- NA;
    }
  }

  output_df[N, ] <- c(gene_cluster,
                      most_frequent_assignment,
                      df$pct_alignment[1],
                      df$cov_in_003[1],
                      df$cov_in_007[1],
                      df$cov_in_011[1],
                      df$cov_in_012[1])

  N = N + 1
}

# store the output data frame as a TAB delimited file so we can import it back
# into the anvi'o additional item data table
write.table(output_df, "gene_clusters_additional_data.txt", quote=FALSE, sep="\t", na="", row.names=FALSE)

