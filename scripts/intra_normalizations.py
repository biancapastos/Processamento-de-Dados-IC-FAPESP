#!/usr/bin/env python3

import pandas as pd
from bioinfokit.analys import norm
import argparse

version = 0.01
parser = argparse.ArgumentParser(description='Intra normalization methods (cpm, tpm and fpkm) script', add_help=True)
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('--matrix', dest='matrix', metavar='matrix.txt',
                    type=str,
                    help='File with integer or combatseq or ruvg matrices to be normalized',
                    required=True)
parser.add_argument('--gene_length_table', dest='gene_length', metavar='gene_length.txt',
                    type=str,
                    help='File with the gene lengths for each gene',
                    required=False)
parser.add_argument('--intra_method', dest='intra_method', metavar='cpm',
                    type=str,
                    help='Normalization method (cpm, tpm, fpkm, none)',
                    required=False, default='cpm')
parser.add_argument('--output_matrix', dest='output_matrix', metavar='output_matrix.txt',
                    type=str, 
                    help='File with the output matrix',
                    required=False)

args = parser.parse_args()
matrix_file = args.matrix
output_matrix_file = args.output_matrix

matrix = pd.read_csv(matrix_file, sep="\t")
matrix = matrix.set_index('Name')

nm = norm()

if args.intra_method == 'cpm' or args.intra_method == 'tpm' or args.intra_method == 'fpkm' or args.intra_method == 'none':
    # CPM
    if args.intra_method == 'cpm':
        nm.cpm(df=matrix)
        cpm_normalized_matrix = nm.cpm_norm
        cpm_normalized_matrix.to_csv(output_matrix_file, sep='\t')
    # gene length for TPM and FPKM normalizations
    elif (args.intra_method == 'tpm') or (args.intra_method == 'fpkm'):
        if not args.gene_length:
            print('Gene length file required for TPM and FPKM normalizations')
            exit(1)
        gene_length_file = args.gene_length
        genelength_table = pd.read_csv(gene_length_file, sep="\t")
        genelength_table = genelength_table.set_index('Name')
        matrix_genelength = pd.merge(matrix, genelength_table, on="Name")
        # TPM
        if args.intra_method == 'tpm':
            nm.tpm(df=matrix_genelength, gl='Length')
            tpm_normalized_matrix = nm.tpm_norm
            tpm_normalized_matrix.to_csv(output_matrix_file, sep='\t')
        # FPKM
        else:
            nm.rpkm(df=matrix_genelength, gl='Length')
            fpkm_normalized_matrix = nm.rpkm_norm
            fpkm_normalized_matrix.to_csv(output_matrix_file, sep='\t')
else:
    print('Norm method not recognized. Please use cpm, tpm, fpkm or none')
    exit(1)