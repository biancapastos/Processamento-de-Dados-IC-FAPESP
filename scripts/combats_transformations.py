#!/usr/bin/env python3

from inmoose.pycombat import pycombat_seq, pycombat_norm
import pandas as pd
import numpy as np
import argparse

version = 0.01
parser = argparse.ArgumentParser(description='Study batch effects calculation for Combat, ComBat_Seq and RemoveBatchEffect', add_help=True)
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('--verbose', dest='verbose', action='store_true')
parser.add_argument('--matrix', dest='matrix', metavar='matrix.txt',
                    type=str,
                    help='File with integer_matrix for combat_seq or any normalizations results matrix for combat and removebatcheffect',
                    required=True)
parser.add_argument('--batch_table', dest='batch_table', metavar='Batch_table.txt',
                    type=str,
                    help='File with the batch effects information for normalization',
                    required=True)
parser.add_argument('--method', dest='method', metavar='method',
                    type=str,
                    help='combat or combatseq or removebatcheffect method',
                    required=True)

parser.add_argument('--output_matrix', dest='output_matrix', metavar='output_matrix.txt',
                    type=str,
                    help='File with combat or combatseq normalizated matrix',
                    required=True)
parser.add_argument('--batch_index', dest='batch_index', metavar='batch_index.txt',
                    type=str,
                    help='File with the study batch effects calculated for removebatcheffect',
                    required=False)

args = parser.parse_args()
matrix_file = args.matrix
batch_table_file = args.batch_table
method = args.method

output_matrix_file = args.output_matrix
studybatch_index_file = args.batch_index

if method == 'combat' or method == 'combatseq' or method == 'removebatcheffect':

    matrix = pd.read_csv(matrix_file, sep="\t")
    matrix = matrix.set_index('Name')

    # if combat method then transform matrix values to arcsinh, else uses integer matrix to combat_seq
    if method == 'combat':
        matrix = matrix.applymap(lambda x: np.arcsinh(x))
        matrix.replace(0, 0.000001, inplace = True)

    if method == 'removebatcheffect':
        matrix.replace(0, 0.000001, inplace = True)
        matrix = matrix.applymap(lambda x: np.log2(x))

    # read table with samples and its different studies(pmid) and tissues
    batch_table = pd.read_table(batch_table_file, sep="\t")

    # calculate in an array the different studies to use as batch effects
    study_batches = batch_table['study'].unique()

    # iterate over the different studies to get the samples from each study in a list and form new dataframes from each study
    list_holder = {}
    df_holder = {}
    for i in range(len(study_batches)):
        list_holder['study' + str(i)] = batch_table.loc[batch_table['study'] == study_batches[i], 'sample'].values.tolist()
        df_holder['study' + str(i)] = pd.DataFrame()

    # populate the dataframes with all genes from the samples from each study
    for i in list_holder:
        df_holder[i] = matrix.loc[:,list_holder[i]]

    # concatenate the matrices in dict to a single dataframe to use as parameter to combat/combat_seq
    study_matrix = pd.DataFrame()
    for i in df_holder.keys():
        aux_matrix = pd.DataFrame(df_holder[str(i)])
        study_matrix = pd.concat([study_matrix,aux_matrix],axis=1)

    # get the batch effects indices to use as parameter for combat/combat_seq function
    batch = []
    datasets = list(df_holder.values())
    for i in range(len(datasets)):
        batch.extend([i for _ in range(len(datasets[i].columns))])

    if method == 'combatseq':
        combatseq_corrected_matrix = pycombat_seq(study_matrix, batch)
        combatseq_corrected_matrix.to_csv(output_matrix_file, sep="\t")
    elif method == 'combat':
        combat_corrected_matrix = pycombat_norm(study_matrix, batch)
        combat_corrected_matrix.to_csv(output_matrix_file, sep="\t")
    elif method == 'removebatcheffect':
        study_matrix.to_csv(output_matrix_file, sep="\t")
        batch_file = pd.DataFrame(batch, columns=['batch'])
        batch_file.to_csv(studybatch_index_file, sep="\t")
else:
    print('Norm method not recognized. Please use combat, combatseq or removebatcheffect')