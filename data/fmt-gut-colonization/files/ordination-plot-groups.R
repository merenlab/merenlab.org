#!/usr/bin/env Rscript

# save this on your disk as 'visualize.R', and then run `chmod +x visualize.R`
# you will run it like this:
#
#    ./visualize.R taxonomy.txt metadata.txt -m mapping_variable -o output_file_prefix --title "title";
#
# for the details of taxonomy.txt, see lab-fiesta #273

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

# command line options
option_list <- list(
		make_option(c("-o", "--output_file_prefix"), default="unknown",
				help = "Output file prefix [default \"%default\"]"),
		make_option(c("-d", "--distance"), default="horn",
				help = "Distance metric [default \"%default\"]"),
		make_option(c("-m", "--mapping_variable"),
				help = "Column in the metadata for sample mapping"),
		make_option("--title", default="(unknown title)",
				help="Title for the output figure [default '%default']")
)

parser <- OptionParser(usage = "script.R [options] input_matrix metadata", option_list=option_list,
		description="A script to generate MDS plots with sample mapping")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

# check if the positional argument is set
if(length(arguments$args) != 2) {
	cat("Incorrect number of required positional arguments\n\n")
	print_help(parser)
	stop()
} else {
	input_file_path <- arguments$args[1]
	metadata_path <- arguments$args[2]
}

if(file.access(input_file_path) == -1){
	stop(sprintf("Input file '%s' does not exist", input_file_path))
}

if(file.access(metadata_path) == -1){
	stop(sprintf("Metadata file '%s' does not exist", metadata_path))
}

if(invalid(options$mapping_variable))
	stop(sprintf("You must define a mapping variable (-m)"))

data <- as.data.frame(read.table(input_file_path, header = TRUE, sep="\t", comment.char = '&'))
metadata <- as.data.frame(read.table(metadata_path, header=TRUE, sep="\t", comment.char = '&'))

if(names(data)[1] != 'key')
	stop(sprintf("Data file '%s' does not seem to be formatted properly", input_file_path))
if(names(metadata)[1] != 'key')
	stop(sprintf("Metadata file '%s' does not seem to be formatted properly", metadata_path))

key_in_both <- intersect(data$key, metadata$key)
data <- data[data$key %in% key_in_both, ]
metadata <- metadata[metadata$key %in% key_in_both, ]

if(dim(metadata)[1] == 0)
	stop("Samples in the input matrix and metadata do not seem to correspond")

metadata <- metadata[match(data$key, metadata$key),]
data$key <- factor(data$key)
metadata$key <- factor(metadata$key)

if(!options$mapping_variable %in% names(metadata)){
	stop(sprintf("Metadata file does not contain mapping variable '%s'", options$mapping_variable))
}

mds <- metaMDS(data[,-1], distance=options$distance)

NMDS = data.frame(MDS1 = mds$points[,1],
                  MDS2 = mds$points[,2],
                  condition=with(metadata, get("condition")),
                  recipient=with(metadata, get("recipient")),
                  first=with(metadata, get("first")),
                  group=with(metadata, get("timing")))

NMDS.mean=aggregate(NMDS[,1:2], list(group=with(metadata, get(options$mapping_variable))), mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}


# data to draw ellipses around groups
df_ell <- data.frame()
only_two_groups <- list()

for(g in levels(NMDS$group)){
    if (nrow(NMDS[NMDS$group==g,]) < 3)
        only_two_groups <- c(only_two_groups, g)
    else
        df_ell <- rbind(df_ell,
                        cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                              veganCovEllipse(cov.wt(cbind(MDS1,MDS2),
                                              wt=rep(1/length(MDS1),
                                              length(MDS2)))$cov,
                                              center=c(mean(MDS1),mean(MDS2))))),
                              group=g))
}

DRAW_SPLIT <- function(){
    if (length(only_two_groups) > 0){
        p <- ggplot(data = NMDS, aes(MDS1, MDS2))
        p <- p + geom_line(aes(x=MDS1, y=MDS2, group=recipient))
		p <- p + geom_point(aes(color = group))
		p <- p + geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=0.5, linetype=1)
		p <- p + geom_line(data = NMDS[NMDS$group %in% only_two_groups,], aes(color = group))
		p <- p + annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group, size=8)
		p <- p + ggtitle(options$title)
		p <- p + theme_bw()
	}
    else{
		p <- ggplot(data = NMDS, aes(MDS1, MDS2))
		p <- p + geom_point(aes(color = group), size=2)
		p <- p + geom_path(data=df_ell, aes(x=MDS1, y=MDS2, colour=group), size=0.5, linetype=1)
		p <- p + geom_polygon(data=df_ell, aes(x=MDS1, y=MDS2, fill=group), alpha=0.1)
		p <- p + facet_grid(recipient ~ .)
		p <- p + ggtitle(options$title)
		p <- p + theme_bw()
	}


    pdf_output <- paste(options$output_file_prefix, "-split.pdf", sep="")
    pdf(pdf_output, width = 6, height = 28, family='Helvetica')
	print(p)
    cat(sprintf("\n\n* PDF for samples split: '%s'\n", pdf_output))
    done <- dev.off()
}


DRAW_TOGETHER <- function(){
    if (length(only_two_groups) > 0){
        p <- ggplot(data = NMDS, aes(MDS1, MDS2))
        p <- p + geom_line(aes(x=MDS1, y=MDS2, group=recipient))
		p <- p + geom_point(aes(color = group))
		p <- p + geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=0.5, linetype=1)
		p <- p + geom_line(data = NMDS[NMDS$group %in% only_two_groups,], aes(color = group))
		p <- p + annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group, size=8)
		p <- p + ggtitle(options$title)
		p <- p + theme_bw()
	}
    else{
		p <- ggplot(data = NMDS, aes(MDS1, MDS2))
		p <- p + geom_line(aes(x=MDS1, y=MDS2, group=recipient), alpha=0.1)
		p <- p + geom_point(aes(color = group), size=4)
		p <- p + geom_path(data=df_ell, aes(x=MDS1, y=MDS2, colour=group), size=0.5, linetype=1)
		p <- p + geom_polygon(data=df_ell, aes(x=MDS1, y=MDS2, fill=group), alpha=0.1)
		p <- p + annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group, size=6)
		p <- p + annotate("text", x=NMDS$MDS1, y=NMDS$MDS2, label=NMDS$first, size=4)
		p <- p + ggtitle(options$title)
		p <- p + theme_bw()

    }


    pdf_output <- paste(options$output_file_prefix, "-together.pdf",sep="")
    pdf(pdf_output, width = 10, height = 6, family='Helvetica')
	print(p)
    cat(sprintf("\n\n* PDF for samples together: '%s'\n\n", pdf_output))
    done <- dev.off()
}

BARS <- function() {
    row.names <- data$key
    col.names <- colnames(data)

    df <- melt(data ,  id = 'key', variable_name = 'bins')

    for (key in levels(df$key)){
      df[df$key == key, ]$value = df[df$key == key, ]$value / sum(df[df$key == key, ]$value) * 100.0
    }

    # only keep taxa that are more than 5% abundant in at least one sample:
    taxa_to_keep <- c();
    i <- 1
    for (taxon in levels(df$bin)) {
      if (max(df[df$bins == taxon, ]$value) > 5)
        taxa_to_keep <- c(taxa_to_keep, taxon)
      i <- i + 1
    }
    df <- df[df$bins %in% taxa_to_keep, ]

    p <- ggplot(df, aes(x=factor(key), y=value, fill=factor(bins)))
    p <- p + geom_bar(position="fill", stat = "identity", width=0.90, colour = 'black')
    p <- p + scale_fill_brewer(palette="Set1")
            # ^^^^ this could be improved to get better colors
            # here is a link for this:
            #     https://www.r-bloggers.com/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
    p <- p + ggtitle(options$title)

    pdf_output <- paste(options$output_file_prefix, "-bars.pdf",sep="")
    pdf(pdf_output, width = 12, height = 6, family='Helvetica')
	print(p)
    cat(sprintf("\n\n* PDF for samples bars: '%s'\n\n", pdf_output))
    done <- dev.off()
}

#####################################################################

DRAW_SPLIT()
DRAW_TOGETHER()
BARS()
