#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:14:00 2023

@author: guohan
"""

import sys
path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Molecule_Analysis'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

import pandas as pd
from similarity_search import main


def run_multiple_similarity_search(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref):
    """

    """
    # read reference SMILES
    df_ref = pd.read_csv(input_file_ref)
    IDs = df_ref[id_column_name_ref].values.tolist()
    SMILES = df_ref[smiles_column_name_ref].values.tolist()
    print(IDs)
    print(SMILES)

    # run similarity search for each reference cmps
    for i, id in enumerate(IDs):
        smiles_column_name = 'Cleaned_SMILES'
        smiles_ref = SMILES[i]
        print(f'Start: {id} {smiles_ref}')
        method = 'substructure'
        output_file = f'{id}.csv'
        main(input_file_library, smiles_column_name, smiles_ref, method, similarity_cutoff=0.5, output_file=output_file, output_option='selected')


if __name__ == '__main__':
    input_file_library = 'library.csv'
    input_file_ref = 'reference_cmps.csv'
    id_column_name_ref = 'ID'
    smiles_column_name_ref = 'Cleaned_Scaffold'
    run_multiple_similarity_search(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref)