#!/usr/bin/env Rscript
library(ggplot2)
library(MASS)
library(ggridges)
library(gridExtra)
library(reshape2)
library(reshape)

# modules that we have identified to be statistically enriched in high-fitness
# donor populations based on recipient colonization patterns after FMT and
# global prevalence.
modules_of_interest <- c("M00570", "M00019", "M00026", "M00432", "M00526", "M00527",
                         "M00018", "M00023", "M00022", "M00015", "M00844", "M00028",
                         "M00048", "M00049", "M00050", "M00051", "M00140", "M00126",
                         "M00924", "M00122", "M00125", "M00115", "M00120", "M00082",
                         "M00083", "M00157", "M00005", "M00007", "M00631", "M00061",
                         "M00855", "M00854", "M00096")

# read the metabolic completion values across genomes
module_completion <- read.table(file='metabolism-completeness-MATRIX.txt', header = TRUE, sep = "\t")
external_genomes <- read.table(file='external-genomes.txt', header = TRUE, sep = "\t")

# let's quickly check the number of genomes per group 
table(external_genomes$cohort)

# turn the boring matrix format into a data frame
df <- melt(module_completion)

# set some meaningful column names
colnames(df) <- c('module', 'genome', 'completion')

# use the external genomes file to associate each genome with a 'group',
# and an individual:
df$group <- external_genomes$cohort[match(df$genome, external_genomes$name)]
df$individual <- external_genomes$individual[match(df$genome, external_genomes$name)]

# subset to the modules of interest
df <- df[df$module %in% modules_of_interest, ]

# some boring steps of defining explicit orders for x-axes in boxplots
HEALTHY_subset <- aggregate(df[df$group == "HEALTHY", 3], list(df[df$group == "HEALTHY", ]$individual), median)
HEALTHY_order <- HEALTHY_subset[order(-HEALTHY_subset$x),]$Group.1
POUCHITIS_subset <- aggregate(df[df$group == "POUCHITIS", 3], list( df[df$group == "POUCHITIS", ]$individual), median)
POUCHITIS_order <- POUCHITIS_subset[order(-POUCHITIS_subset$x),]$Group.1
CROHNS_subset <- aggregate(df[df$group == "CROHNS", 3], list( df[df$group == "CROHNS", ]$individual), median)
CROHNS_order <- CROHNS_subset[order(-CROHNS_subset$x), ]$Group.1
individuals_order <- c(c('FMT_HIGH_FITNESS', 'FMT_LOW_FITNESS'), HEALTHY_order, POUCHITIS_order, CROHNS_order)

# set explicit group orders, and assign some group colors
groups_order <- c("FMT_HIGH_FITNESS", "FMT_LOW_FITNESS", "HEALTHY", "POUCHITIS", "CROHNS")
group_colors <- c("#ec5f1c", "#034f84", "#feb236", "#86af49", "#ff0202")

# plot the boxplots
pdf(file = "boxplots.pdf",  width = 13, height = 5)
plot_individuals <- ggplot(data=df, aes(x=individual, y=completion, group=individual)) +
  geom_boxplot(aes(fill=group), alpha=0.35, outlier.shape = NA, color='#808080') +
  geom_jitter(colour='#222222', width = 0.20, height = 0.02, size=0.1, alpha=0.5) +
  theme_bw() +
  theme(legend.position="bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("Individuals") +
  scale_x_discrete(limits = individuals_order) +
  scale_fill_manual(values = group_colors)

plot_groups <- ggplot(data=df, aes(x=group, y=completion, group=group)) +
  geom_boxplot(aes(fill=group), alpha=0.35, outlier.shape = NA, color=NA) +
  geom_violin(fill="#505050", alpha=0.35, width=1.3, colour = '#505050') +
  geom_jitter(colour='#222222', width = 0.3, height = 0.02, size=0.1, alpha=0.05) +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Metabolic module completion") +
  xlab("Groups") +
  scale_x_discrete(limits = groups_order) +
  scale_fill_manual(values = group_colors)

grid.arrange(plot_groups, plot_individuals, ncol=2, widths=c(1, 3))
dev.off()

# some of the most significant modules that differ between healthy group vs the
# IBD group to show some examples
modules_to_display <- c("M00924","M00122", "M00023", "M00028", "M00570", "M00082", "M00844", "M00015",  "M00526", "M00022") 
dfx <- df[df$module %in% modules_to_display, ]
dfx$module = factor(dfx$module, levels=modules_to_display)


# some boring steps of defining explicit orders for x-axes in boxplots
HEALTHY_subsetx <- aggregate(dfx[dfx$group == "HEALTHY" & dfx$module == "M00924", 3], list(dfx[dfx$group == "HEALTHY" & dfx$module == "M00924", ]$individual), mean)
HEALTHY_orderx <- HEALTHY_subsetx[order(-HEALTHY_subsetx$x),]$Group.1
POUCHITIS_subsetx <- aggregate(dfx[dfx$group == "POUCHITIS" & dfx$module == "M00924", 3], list( dfx[dfx$group == "POUCHITIS" & dfx$module == "M00924", ]$individual), mean)
POUCHITIS_orderx <- POUCHITIS_subsetx[order(-POUCHITIS_subsetx$x),]$Group.1
CROHNS_subsetx <- aggregate(dfx[dfx$group == "CROHNS" & dfx$module == "M00924", 3], list( dfx[dfx$group == "CROHNS" & dfx$module == "M00924", ]$individual), mean)
CROHNS_orderx <- CROHNS_subsetx[order(-CROHNS_subsetx$x), ]$Group.1
individuals_orderx <- c(c('FMT_HIGH_FITNESS', 'FMT_LOW_FITNESS'), HEALTHY_orderx, POUCHITIS_orderx, CROHNS_orderx)


pdf(file = "ridges-for-metabolisms.pdf",  width = 13, height = 5)
ggplot(data=dfx, aes(x = completion, y = individual, fill = group)) +
  geom_density_ridges(scale = 10, size = 0.25, rel_min_height = 0.01, alpha=0.35, colour="#222222") +
  theme_ridges() +
  scale_y_discrete(limits = individuals_orderx) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(breaks=c(0, 1)) +
  coord_cartesian(xlim = c(-0.5, 1.5)) +
  theme(legend.position="bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_wrap(. ~ module, nrow = 2)
dev.off()


pdf(file = "boxplots-for-metabolisms.pdf",  width = 13, height = 5)
ggplot(data=dfx, aes(x = completion, y = group, fill = group)) +
  geom_boxplot(aes(fill=group), alpha=0.35, outlier.shape = NA, color='#505050') +
  geom_jitter(colour='#222222', width = 0.1, height = 0.1, size=0.5, alpha=0.1) +
  theme_ridges() +
  scale_y_discrete(limits = c('FMT_HIGH_FITNESS', 'FMT_LOW_FITNESS', 'HEALTHY', 'POUCHITIS', 'CROHNS')) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(breaks=c(0, 1)) +
  coord_cartesian(xlim = c(-0, 1)) +
  theme(legend.position="bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_wrap(. ~ module, nrow = 2)
dev.off()


# print some stats for ridges
printf <- function(...) invisible(print(sprintf(...)))
printf("DIFFERENCES BETWEEN HEALTHY vs CROHNS + POUCHITIS for RIDGELINE PLOTS")
for(module in modules_to_display){
    w <- wilcox.test(dfx[dfx$group == "HEALTHY" & dfx$module == module, ]$completion, dfx[dfx$group %in% c("CROHNS", "POUCHITIS") & dfx$module == module, ]$completion)
    printf("%s: %f", module, w$p.value)
}


## testing whether differences are significant
wilcox.test(df[df$group == "FMT_HIGH_FITNESS", ]$completion, df[df$group == "FMT_LOW_FITNESS", ]$completion)
wilcox.test(df[df$group == "FMT_HIGH_FITNESS", ]$completion, df[df$group == "POUCHITIS", ]$completion)
wilcox.test(df[df$group == "FMT_HIGH_FITNESS", ]$completion, df[df$group == "CROHNS", ]$completion)
wilcox.test(df[df$group == "POUCHITIS", ]$completion, df[df$group == "CROHNS", ]$completion)
wilcox.test(df[df$group == "FMT_HIGH_FITNESS", ]$completion, df[df$group == "HEALTHY", ]$completion)

# mean genome completion
mean(external_genomes[external_genomes$cohort == "HEALTHY", ]$genome_completion)
mean(external_genomes[external_genomes$cohort == "POUCHITIS", ]$genome_completion)
mean(external_genomes[external_genomes$cohort == "CROHNS", ]$genome_completion)
mean(external_genomes[external_genomes$cohort == "FMT_LOW_FITNESS", ]$genome_completion)
mean(external_genomes[external_genomes$cohort == "FMT_HIGH_FITNESS", ]$genome_completion)

# mean genome length
mean(external_genomes[external_genomes$cohort == "HEALTHY", ]$total.length)
mean(external_genomes[external_genomes$cohort == "POUCHITIS", ]$total.length)
mean(external_genomes[external_genomes$cohort == "CROHNS", ]$total.length)
mean(external_genomes[external_genomes$cohort == "CROHNS" | external_genomes$cohort == "POUCHITIS", ]$total.length)
mean(external_genomes[external_genomes$cohort == "FMT_LOW_FITNESS", ]$total.length)
mean(external_genomes[external_genomes$cohort == "FMT_HIGH_FITNESS", ]$total.length)

# num invdividuals per group
length(unique(external_genomes[external_genomes$cohort == "HEALTHY", ]$individual))
length(unique(external_genomes[external_genomes$cohort == "POUCHITIS", ]$individual))
length(unique(external_genomes[external_genomes$cohort == "CROHNS", ]$individual))

mean(table(external_genomes[external_genomes$cohort == "HEALTHY", ]$individual))
mean(table(external_genomes[external_genomes$cohort == "POUCHITIS", ]$individual))
mean(table(external_genomes[external_genomes$cohort == "CROHNS", ]$individual))
mean(table(external_genomes[external_genomes$cohort == "FMT_LOW_FITNESS", ]$individual))
mean(table(external_genomes[external_genomes$cohort == "FMT_HIGH_FITNESS", ]$individual))

