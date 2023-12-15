#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:14:00 2023

@author: guohan
"""

import sys, time, os
path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Molecule_Analysis'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_size(10)
from similarity_search import main
from utils.tools import remove_unnamed_columns


def run_multiple_similarity_search(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref,
                                   method, similarity_cutoff, analogs_dir):
    """
    Perform similarity search for multiple SMILES from 'input_file_ref'
    """
    # read reference SMILES
    df_ref = pd.read_csv(input_file_ref)
    IDs = df_ref[id_column_name_ref].values.tolist()
    SMILES = df_ref[smiles_column_name_ref].values.tolist()
    print(IDs)
    print(SMILES)

    # run similarity search for each reference cmps
    total_time = []
    for i, id in enumerate(IDs):
        start_time = time.time()
        smiles_column_name = 'Cleaned_SMILES'
        smiles_ref = SMILES[i]
        print(f'Start: {id} {smiles_ref}')
        output_folder = analogs_dir
        output_file = f'{id}.csv'
        main(input_file_library, smiles_column_name, smiles_ref, method, similarity_cutoff=similarity_cutoff, output_folder=output_folder, output_file=output_file, output_option='satisfied')
        end_time = time.time()
        elapsed_time = end_time - start_time
        total_time.append(elapsed_time)
        print(f'Done, time = {elapsed_time:.2f}')
    print(f'Mean time is {np.mean(np.array(total_time))}, std is {np.std(np.array(total_time))}')


def get_best_score(input_file, analogs_dir, n=0):
    """
    Get the top-n most similar compound
    :param input_file: str, path of the file containing reference compounds
    :param analogs_dir: str, path of the folder containing analogs
    :param n: int, indicate the top-n most similar analog
    """
    # output file
    output_file = os.path.splitext(os.path.abspath(input_file))[0]

    df = pd.read_csv(input_file)
    df_bestScore = pd.DataFrame(df, columns = ['ID', 'SMILES', 'Cleaned_SMILES'])

    df_bestScore['Analog'] = df_bestScore['ID'].apply(lambda cmp: get_analog(cmp, analogs_dir, n))
    df_bestScore[['Analog_ID', 'Analog_SMILES', 'Analog_Cleaned_SMILES', 'Similarity_Score']] = \
        pd.DataFrame(df_bestScore['Analog'].tolist())
    COLUMNS = ['ID', 'SMILES', 'Cleaned_SMILES', 'Analog_ID', 'Analog_SMILES', 'Analog_Cleaned_SMILES', 'Similarity_Score']
    df_bestScore = pd.DataFrame(df_bestScore, columns = COLUMNS)

    df_bestScore = df_bestScore.reset_index(drop=True)
    print('Number of rows:', df_bestScore.shape[0])
    df_bestScore = remove_unnamed_columns(df_bestScore)
    df_bestScore.to_csv(f'{output_file}_Top{n+1}_{df_bestScore.shape[0]}.csv')


def get_analog(cmp, analogs_dir, n):
    files = os.listdir(analogs_dir)
    filtered_files = [file for file in files if file.startswith(f'{cmp}_')]
    df_cmp = pd.read_csv(f'{analogs_dir}/{filtered_files[0]}')

    if df_cmp.shape[0] <= n:
        return ['', '', '', 0.0]

    return df_cmp.iloc[n, [False, True, True, True, False, False, False, False, False, True]].tolist()


def plot_distribution(input_file):
    """

    """
    # output file
    output_file = os.path.splitext(os.path.abspath(input_file))[0]

    df = pd.read_csv(input_file)

    # plot distribution
    values = df['Similarity_Score'].tolist()
    output_file = f'{output_file}_similarity_distribution.pdf'

    plt.figure(1)
    plt.hist(values, 10, range=[0.0, 1.0])
    plt.xlabel('Similarity Score', fontproperties=font)
    plt.ylabel('Counts', fontproperties=font)
    plt.xticks(fontproperties=font)
    plt.yticks(fontproperties=font)

    plt.savefig(output_file, dpi=300)
    plt.show()
    plt.close()


if __name__ == '__main__':
    input_file_library = 'library.csv'
    # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'

    input_file_ref = 'reference_cmps.csv'
    id_column_name_ref = 'ID'
    smiles_column_name_ref = 'Cleaned_SMILES'

    method = 'mcs'
    similarity_cutoff = 0.3
    analogs_dir = 'library_FP0.3'
    run_multiple_similarity_search(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref,
                                   method, similarity_cutoff, analogs_dir)