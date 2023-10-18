#!/usr/bin/env python3

import pandas as pd
import argparse

version = 0.01
parser = argparse.ArgumentParser(description='Rounds estimated counts from Salmon to integer and removal of non-expressed genes', add_help=True)
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('--verbose', dest='verbose', action='store_true')
parser.add_argument('--salmon_numreads', dest='salmon_numreads', metavar='salmon_numreads_matrix.txt',
                    type=str,
                    help='File with the estimated counts from Salmon quantmerge',
                    required=True)
parser.add_argument('--integer_matrix', dest='integer_matrix', metavar='integer_matrix.txt',
                    type=str,
                    help='File with the Salmon numreads matrix to integer values',
                    required=True)

args = parser.parse_args()
salmon_numreads_file = args.salmon_numreads
integer_matrix_file = args.integer_matrix

numreads_matrix = pd.read_csv(salmon_numreads_file, sep="\t")
numreads_matrix = numreads_matrix.set_index('Name')
numreads_rounded_matrix = numreads_matrix.round(0).astype(int)
integer_matrix_wo_zeros = numreads_rounded_matrix.drop(numreads_rounded_matrix[(numreads_rounded_matrix.iloc[:, :] == 0).all(axis=1)].index)
integer_matrix_wo_zeros.to_csv(integer_matrix_file, sep="\t")