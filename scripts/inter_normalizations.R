#!/usr/bin/env Rscript

library(argparse)
library(edgeR)
library(preprocessCore)

parser <- ArgumentParser(description= 'Inter Normalization methods (tmm, ctf, uq, cuf)')

parser$add_argument('--matrix', dest='matrix', 
                    help= 'matrix to be normalized')
parser$add_argument('--norm_method', dest='norm_method', default= 'tmm',
                    help= 'inter normalization method (tmm, ctf, uq, cuf)')
parser$add_argument('--output_matrix', dest='output_matrix', 
                    help= 'output normalized matrix ')

xargs<- parser$parse_args()

if(xargs$norm_method == 'tmm' || xargs$norm_method == 'ctf' || xargs$norm_method == 'uq' || xargs$norm_method == 'cuf'){

    x <- read.delim(xargs$matrix, row.names = "Name", header=TRUE)
    y = DGEList(counts = x)

    sink(xargs$output_matrix)
    cat("Name")
    sink()

    if( xargs$norm_method == 'tmm' || xargs$norm_method == 'ctf' ){
        
        norm_factors <- calcNormFactors(y, method="TMM")

        if(xargs$norm_method == 'tmm'){
            TMM_normalized <- sweep(y$counts, 2, c((norm_factors$samples$lib.size*norm_factors$samples$norm.factors)/(10**6)), "/")
            TMM_rounded <- apply(TMM_normalized, 2, round, digits = 5)
            as.data.frame(TMM_rounded)
            write.table(TMM_rounded, xargs$output_matrix, append=TRUE, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
        }else if(xargs$norm_method == 'ctf'){
            CTF_normalized <- sweep(y$counts, 2, c(norm_factors$samples$norm.factors), "/")
            CTF_rounded <- apply(CTF_normalized, 2, round, digits = 5)
            as.data.frame(CTF_rounded)
            write.table(CTF_rounded, xargs$output_matrix, append=TRUE, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
        }

    }else if(xargs$norm_method == 'uq' || xargs$norm_method == 'cuf'){

        norm_factors <- calcNormFactors(y, method="upperquartile")

        if(xargs$norm_method == 'uq'){
            UQ_normalized <- sweep(y$counts, 2, c((colSums(y$counts)*norm_factors$samples$norm.factors)/(10**6)), "/")
            UQ_rounded <- apply(UQ_normalized, 2, round, digits = 5)
            as.data.frame(UQ_rounded)
            write.table(UQ_rounded, xargs$output_matrix, append=TRUE, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
        }else if(xargs$norm_method == 'cuf'){
            CUF_normalized <- sweep(y$counts, 2, c(norm_factors$samples$norm.factors), "/")
            CUF_rounded <- apply(CUF_normalized, 2, round, digits = 5)
            as.data.frame(CUF_rounded)
            write.table(CUF_rounded, xargs$output_matrix, append=TRUE, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
        }
    }

}else{
    print('Norm method not recognized. Please use tmm, ctf, uq or cuf')
}

