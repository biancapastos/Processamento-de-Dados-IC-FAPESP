#!/usr/bin/env Rscript

library(preprocessCore)
library(argparse)
library(RUVSeq)

parser <- ArgumentParser(description= 'Removal of Unwanted Variation RUV normalization method - RUVg approach')

parser$add_argument('--integer_matrix', dest='integer_matrix', 
                    help= 'matrix to be normalized')
parser$add_argument('--spikes', dest='spikes', 
                    help= 'spikes (negative control genes) to use as parameter')
parser$add_argument('--groups', dest='groups', 
                    help= 'groups to use as parameter')

parser$add_argument('--output_matrix', dest='output_matrix', 
                    help= 'output normalized matrix ')

xargs<- parser$parse_args()

matrix <- read.csv(xargs$integer_matrix, header=TRUE, sep='\t')
genes <- matrix$Name
row.names(matrix) <- matrix$Name
matrix <- subset(matrix, select = -Name)

spikes <- read.csv(xargs$spikes, header=TRUE, sep='\t')
#spikes <- spikes$Name

groups <- read.csv(xargs$groups, header=TRUE, sep='\t')

# Cluster data into corrected groups
x <- as.factor(groups$group)
set <- newSeqExpressionSet(as.matrix(matrix),
                           phenoData = data.frame(x, row.names=colnames(matrix)))

# Upperquartile normalization
set <- betweenLaneNormalization(set, which="upper")

# RUVg normalization
set1 <- RUVg(set, spikes$Name, k=1)
norm_values <- normCounts(set1)

sink(xargs$output_matrix)
cat("Name")
sink()
write.table(norm_values, xargs$output_matrix, append=TRUE, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
