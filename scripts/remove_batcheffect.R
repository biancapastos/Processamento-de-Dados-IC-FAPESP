#!/usr/bin/env Rscript

library(preprocessCore)
library(argparse)
library(limma)

parser <- ArgumentParser(description= 'RemoveBatchEffect normalization method for log2 matrix')

parser$add_argument('--matrix', dest='matrix', 
                    help= 'matrix to be normalized')
parser$add_argument('--batch_effect', dest='batch_effect', 
                    help= 'batch_effect index to use as parameter')

parser$add_argument('--output_matrix', dest='output_matrix', 
                    help= 'output normalized matrix ')

xargs<- parser$parse_args()

study_matrix <- read.delim(xargs$matrix, row.names = 'Name', header=TRUE)
study_batch_effects <- read.delim(xargs$batch_effect, header=TRUE, sep='\t')
study_batch_effects <- study_batch_effects[,2]
study_matrix_adjusted = removeBatchEffect(study_matrix, batch=study_batch_effects)

sink(xargs$output_matrix)
cat("Name")
sink()
write.table(study_matrix_adjusted, xargs$output_matrix, append=TRUE, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")