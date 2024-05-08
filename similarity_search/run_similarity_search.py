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

from similarity_search import similarity_search_single_query, similarity_search_multiple_query, select_analogs, plot_distribution


def run_similarity_search_single_query():
    """
    Perform similarity search for single query SMILES
    """
    input_file_library = 'tests/library.csv'
    # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
    smiles_column_name = 'Cleaned_SMILES'
    smiles_query = 'N#CC1(Oc2ccc(-c3ccccc3)cc2)CC2(CCN(C(=O)Nc3ccccc3)C2)C1'
    # smiles_query = 'C1=CC=CC=C1'

    output_folder = 'tests/similarity_search_results'
    output_file = 'output.csv'
    output_option = 'satisfied'

    ### Method: 'fingerprint' ###
    # similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method='fingerprint', similarity_cutoff=0.3, output_folder=output_folder, output_file=output_file, output_option=output_option)
    ### Method: 'mcs' ###
    similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method='mcs', mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_file=output_file, output_option=output_option)
    ### Method: 'substructure' ###
    # similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method='substructure', substructure_method='SMARTS', output_folder=output_folder, output_file=output_file, output_option=output_option)


def run_similarity_search_multiple_query():
    """
    Perform similarity search for multiple query SMILES in input_file_query
    """
    input_file_library = 'tests/library.csv'
    # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
    input_file_query = 'tests/query_cmps.csv'
    id_column_name_query = 'ID'
    smiles_column_name_query = 'Cleaned_SMILES'

    output_folder = 'tests/similarity_search_results'
    output_option = 'satisfied'

    ### Method: 'fingerprint' ###
    # similarity_search_multiple_query(input_file_library, input_file_query, id_column_name_query, smiles_column_name_query,
    #                                method='fingerprint', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
    ### Method: 'mcs' ###
    similarity_search_multiple_query(input_file_library, input_file_query, id_column_name_query, smiles_column_name_query,
                                   method='mcs', mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
    ### Method: 'substructure' ###
    # similarity_search_multiple_query(input_file_library, input_file_query, id_column_name_query, smiles_column_name_query,
    #                                method='substructure', substructure_method='SMARTS',output_folder=output_folder, output_option=output_option)


def run_select_analogs():
    """
    Select analogs
    """
    input_file_query = 'tests/query_cmps.csv'
    analogs_dir = 'tests/similarity_search_results'
    ### Select analogs above similarity cutoff
    # analog_method = 'cutoff'
    # similarity_cutoff = 0.4
    # select_analogs(input_file_query, analogs_dir, analog_method, similarity_cutoff=similarity_cutoff, deduplication=False)
    ### Select the rank-n most similar analog
    analog_method = 'rank'
    similarity_rank = 1
    select_analogs(input_file_query, analogs_dir, analog_method, similarity_rank=similarity_rank, deduplication=False)


def run_plot_distribution():
    input_file = 'tests/query_cmps_Top1_12.csv'
    plot_distribution(input_file)



if __name__ == '__main__':
    # run_similarity_search_single_query()
    run_similarity_search_multiple_query()
    # run_select_analogs()
    # run_plot_distribution()