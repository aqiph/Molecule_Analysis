#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 10:12:33 2022

@author: guohan

1. calculate the number of scaffolds in the given input file
2. plot distance distribution based on fingerprint

"""

import os
import pandas as pd
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_size(12)

from utils.molecular_description import get_fingerprint, cal_fingerprint_distance, get_scaffold



### calculate the number of scaffolds in the given input file ###

def cal_num_scaffold(input_file, smiles_column_name = 'SMILES'):
    """
    calculate the number of scaffolds in the given library
    :param input_file: str, the filename of the input file
    :param smiles_column_name: str, the name of the SMILES column
    """
    # output file path without extension
    output_file, fmt = os.path.splitext(os.path.abspath(input_file))

    if fmt in {'.csv'}:
        df = pd.read_csv(input_file)
    elif fmt in {'.xlsx'}:
        df = pd.read_excel(input_file)
    else:
        print('Error: Invalid format of the input file')
        return

    # get scaffolds in the given input file
    fn = partial(get_scaffold, include_chirality = False, generic = False)
    df['Scaffold'] = df[smiles_column_name].apply(fn)
    df_scaffold = pd.DataFrame(df, columns=['Scaffold'])

    # remove duplicates scaffolds
    df_scaffold.drop_duplicates(subset = 'Scaffold', keep='first', inplace=True, ignore_index=True)
    df_scaffold['ID'] = np.arange(0, df_scaffold.shape[0], 1)
    df_scaffold = df_scaffold[['ID', 'Scaffold']]

    # write to file
    df_scaffold = df_scaffold.reset_index(drop=True)
    print('Number scaffolds:', df_scaffold.shape[0])
    df_scaffold.to_csv('{}_scaffold.csv'.format(output_file))


### plot distance distribution based on fingerprint  ###

def plot_distance_distribution_by_fingerprint(input_file, smiles_column_name = 'SMILES', fp_method = 'ecfp4'):
    """
    calculate the distance between every pair of SMILES, based on fingerprint
    plot the distribution of the distance
    :param input_file: str, the filename of the input file
    :param smiles_column_name: str, the name of the SMILES column
    :param fp_method: str, method to compute fingerprint (topology: topological fingerprint;
                                                         maccs: MACCS fingerprint;
                                                         atomPairs: atom pairs;
                                                         torsions: topological torsions;
                                                         ecfp4: Morgan ECFP4 (connectivity) fingerprint;
                                                         fcfp4: Morgan FCFP4 (feature-based) fingerprint)
                                                         mhfp6: MinHash fingerprint
    """
    # output file path without extension
    output_file, fmt = os.path.splitext(os.path.abspath(input_file))

    if fmt in {'.csv'}:
        df = pd.read_csv(input_file)
    elif fmt in {'.xlsx'}:
        df = pd.read_excel(input_file)
    else:
        print('Error: Invalid format of the input file')
        return

    # calculate fingerprint
    fp_method = fp_method.lower()

    smiles_list = df[smiles_column_name].values.tolist()
    num_smiles = len(smiles_list)
    fp_list = [get_fingerprint(smiles, fp_method) for smiles in smiles_list]

    # calculate distance between every pair of SMILES
    distances = []

    for i in range(num_smiles):
        for j in range(i):
            d = cal_fingerprint_distance(fp_list[i], fp_list[j], fp_method)
            distances.append(d)

    if True:
        print('maximum:', max(distances), 'minimum:', min(distances))

    # plot distance distribution
    output_file_distribution = '{}_distance_distribution.pdf'.format(output_file)

    plt.figure(1)
    plt.hist(distances, 50, range=[0, 1.0])
    plt.xlabel('Fingerprint Distance', fontproperties=font)
    plt.ylabel('Counts', fontproperties=font)
    plt.xticks(fontproperties=font)
    plt.yticks(fontproperties=font)

    plt.savefig(output_file_distribution, dpi=300)
    plt.show()
    plt.close()



if __name__ == '__main__':

    ### calculate the number of scaffolds in the given input file ###
    input_file = 'diversity_analysis/tests/test_diversity_analysis.csv'
    smiles_column_name = 'Cleaned_SMILES'
    cal_num_scaffold(input_file, smiles_column_name)

    ### calculate the number of scaffolds in the given input file ###
    input_file = 'diversity_analysis/tests/test_diversity_analysis.csv'
    smiles_column_name = 'Cleaned_SMILES'
    plot_distance_distribution_by_fingerprint(input_file, smiles_column_name, fp_method='ecfp4')
