#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)

#docker run --rm -it -v $PWD:/input rocker/tidyverse:4.3.1 /bin/bash


##########################################################
## inputs
###########################################################
args <- commandArgs(trailingOnly = TRUE)
fileout <- args[length(args)] # Letztes Argument ist die Ausgabedatei
filein <- args[-length(args)] # Alle anderen Argumente sind Eingabedateien







# * Alignment Coverage Plot:
#     * Show the alignment start and end positions of the query sequence against the subject sequence. This helps visualize how much of the query sequence is covered in each alignment.
# * Percent Identity Distribution:
#     * Create a histogram of percent identity across all matches. This can give insights into how similar the matched sequences are.
# * Alignment Length vs. Percent Identity:
#     * Use a scatter plot to explore the relationship between alignment length and percent identity. Larger alignments with high identity often indicate more significant matches.
# * E-value Distribution:
#     * Plot a histogram or boxplot of e-values to identify the quality of matches. Lower e-values signify more reliable hits.
# * Query Coverage Heatmap:
#     * Build a heatmap showing the regions of the query sequence covered by subject sequences. This is useful for spotting gaps or overlaps.
# * Bit Score vs. Percent Identity:
#     * Scatter plot or bubble plot showing the relationship between bit score and percent identity for matched sequences.
# * Hit Locations:
#     * Visualize where hits occur on the query or subject sequence, possibly mapping overlaps or gaps.