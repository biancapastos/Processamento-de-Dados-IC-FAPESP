#!/usr/bin/env Rscript

library(argparse)
library(edgeR)
library(preprocessCore)

parser <- ArgumentParser(description= 'qnt inter normalization method')

parser$add_argument('--matrix', dest='matrix', 
                    help= 'matrix to be normalized')

parser$add_argument('--output_matrix', dest='output_matrix', 
                    help= 'output normalized matrix ')

xargs<- parser$parse_args()

x <- read.delim(xargs$matrix, row.names = "Name", header=TRUE)
y = DGEList(counts = x)

sink(xargs$output_matrix)
cat("Name")
sink()

normalized_data <- normalize.quantiles(y$counts)
colnames(normalized_data) <- colnames(x)
rownames(normalized_data) <- rownames(x)
QNT_rounded <- apply(normalized_data, 2, round, digits = 5)
as.data.frame(QNT_rounded)
write.table(QNT_rounded, xargs$output_matrix, append=TRUE, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
