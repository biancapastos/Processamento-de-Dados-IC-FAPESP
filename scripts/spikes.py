#!/usr/bin/env python3

import pandas as pd
import numpy as np

from scipy import stats
from scipy.stats import shapiro
import argparse

version = 0.01
parser = argparse.ArgumentParser(description='Calculation of negative control genes (spikes) for ruvg', add_help=True)
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('--verbose', dest='verbose', action='store_true')
parser.add_argument('--tpm_matrix', dest='tpm_matrix', metavar='tpm_matrix.txt',
                    type=str,
                    help='File with the tpm normalized matrix from integer matrix',
                    required=True)

parser.add_argument('--spikes', dest='spikes', metavar='spikes.txt',
                    type=str,
                    help='File with the genes considered spikes',
                    required=True)

args = parser.parse_args()
tpm_matrix_file = args.tpm_matrix
spikes_file = args.spikes

tpm_matrix = pd.read_csv(tpm_matrix_file, sep="\t")
tpm_matrix = tpm_matrix.set_index('Name')

# calculation of CoV, mean and standard deviation (sd)
cov_values = (tpm_matrix.std(axis=1) / tpm_matrix.mean(axis=1)) * 100
mean_expression = tpm_matrix.mean(axis=1)
sd_expression = tpm_matrix.std(axis=1)
# Make new dataframe for above statistics calculated
stats_df = pd.DataFrame({"CoV": cov_values, "Media": mean_expression, "DevPad": sd_expression})
stats_df.index = tpm_matrix.index

# Shapiro test
def shapiro_pvalue(expression):
    shapiro_stat, shapiro_pvalue = stats.shapiro(expression)
    return shapiro_pvalue

stats_df['Shapiro'] = tpm_matrix.apply(shapiro_pvalue, axis=1)

#Index (genes names) go through 3 filters:
# 1) CoV < 25%;
# 2) Samples means > 100 (mean genes expressed);
# 3) Shapiro Wilk's Teste_p < que 0.1.
spikes = stats_df[(stats_df['CoV']<25)&
                      (stats_df['Media']>100)&
                      (stats_df['Shapiro']<0.1)]['Shapiro'].index
spikes = spikes.to_frame(index=False)
spikes.to_csv(spikes_file, sep="\t")