#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:10:00 2023

@author: guohan

1. Similarity search
2. Select analogs based on similarity cutoff or ranking
3. Plot similarity score distribution

"""


import os, time
from tqdm import tqdm
import pandas as pd
import numpy as np
from rdkit import Chem

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_size(10)

from utils.molecular_description import get_fingerprint, cal_fingerprint_distance, cal_MCS
from utils.tools import remove_unnamed_columns



### Similarity search ###
def is_similar_by_fingerprint(smiles_1, smiles_2, similarity_cutoff=0.7):
    """
    Determine whether smiles_1 is similar to smiles_2 based on ecfp4 fingerprint
    :param smiles_1: str, SMILES 1
    :param smiles_2: str, SMILES 2
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles_1 is considered to be similar to smiles_2
    :return: bool, whether smiles_1 is considered to be similar to smiles_2
    """
    fp_1 = get_fingerprint(smiles_1, fp_method = 'ecfp4')
    fp_2 = get_fingerprint(smiles_2, fp_method = 'ecfp4')
    dist = cal_fingerprint_distance(fp_1, fp_2, fp_method='ecfp4')
    similarity_score = 1.0 - dist
    return similarity_score >= similarity_cutoff, similarity_score


def is_similar_by_mcs(smiles_1, smiles_2, match='exact', similarity_cutoff=0.7):
    """
    Determine whether smiles_1 is similar to smiles_2 based on MCS
    :param smiles_1: str, SMILES 1
    :param smiles_2: str, SMILES 2
    :param match: str, specify 'match' to use different comparison functions, allowed values include 'exact', 'anyAtom'
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles_1 is considered to be similar to smiles_2
    :return: bool, whether smiles_1 is considered to be similar to smiles_2
    """
    similarity_score = cal_MCS(smiles_1, smiles_2, match)
    return similarity_score >= similarity_cutoff, similarity_score


def is_similar_by_substructure(smiles, smiles_sub, substructure_method='SMILES'):
    """
    Determine whether smiles_sub is a substructure of smiles
    :param smiles: str, SMILES
    :param smiles_sub: str, SMILES for substructure
    :param substructure_method: str, method used to convert smiles_sub to mol_sub, allowed values include 'SMILES', 'SMARTS'
    :return: bool, whether smiles_sub is a substructure of smiles
    """
    # check SMILES
    if not smiles:
        print(f"Error: Invalid SMILES {smiles}")
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        print(f"Error: Invalid SMILES {smiles}")
        return False
    if not mol:
        print(f"Error: Invalid SMILES {smiles}")
        return False

    # determine whether smiles_sub is a substructure of smiles
    if substructure_method == 'SMILES':
        mol_sub = Chem.MolFromSmiles(smiles_sub)
    elif substructure_method == 'SMARTS':
        smiles_sub = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_sub))   # Change molecule before get mol as a pattern, not using the Kekule form
        # print(smiles_sub)
        mol_sub = Chem.MolFromSmarts(smiles_sub)
    else:
        raise Exception('Error: Invalid substructure method.')

    hasSubstructure = mol.HasSubstructMatch(mol_sub)
    return hasSubstructure, int(hasSubstructure)


def similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method, **kwargs):
    """
    Find SMILES in input_file_library that are similar to the query SMILES.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param input_file_library: str, path of the input library.
    :param smiles_column_name: str, name of the SMILES column.
    :param smiles_query: str, query SMILES.
    :param method: str, method for similarity search, allowed values include 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the query compound.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param mcs_match: str, specify 'mcs_match' to use different comparison functions for MCS, allowed values include 'exact', 'anyAtom'.
    :param substructure_method: str, method used to convert smiles_query to mol_query, allowed values include 'SMILES', 'SMARTS'.
    :param output_folder: str, path of the output file.
    :param output_file: str, name of the output file.
    :param output_option: str, options to output data, allowed values include 'satisfied', 'not_satisfied' and 'all'.
    """
    # files
    df = pd.read_csv(input_file_library)
    df['ID'] = df['ID'].apply(lambda id: str(id))
    COLUMNS = df.columns.tolist() + ['Similarity_Score']
    assert smiles_column_name in COLUMNS, 'Error: Invalid input SMILES column name.'

    folder = kwargs.get('output_folder', os.getcwd())
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f'Fold {folder} created.')
    output_file = kwargs.get('output_file', input_file_library)
    output_file = os.path.splitext(os.path.split(output_file)[1])[0]
    output_file = os.path.join(folder, output_file)

    # check query SMILES
    if not smiles_query:
        print(f"Error: Invalid SMILES {smiles_query}")
        return None
    try:
        mol = Chem.MolFromSmiles(smiles_query)
    except Exception:
        print(f"Error: Invalid SMILES {smiles_query}")
        return None
    if not mol:
        print(f"Error: Invalid SMILES {smiles_query}")
        return None

    # determine whether SMILES in df is similar to the query SMILES
    method = method.lower()
    if method == 'fingerprint':
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_fingerprint(smiles, smiles_query, similarity_cutoff))
    elif method == 'mcs':
        match = kwargs.get('mcs_match', 'exact')
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_mcs(smiles, smiles_query, match, similarity_cutoff))
    elif method == 'substructure':
        substructure_method = kwargs.get('substructure_method', 'SMILES')
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_substructure(smiles, smiles_query, substructure_method))
    else:
        raise Exception('Error: Invalid similarity search method.')
    df[['is_similar', 'Similarity_Score']] = pd.DataFrame(df['is_similar'].values.tolist())

    # prepare output file
    output_option = kwargs.get('output_option', 'satisfied')
    if output_option == 'satisfied':
        df_sim = pd.DataFrame(df[df['is_similar']], columns = COLUMNS)
    elif output_option == 'not_satisfied':
        df_sim = pd.DataFrame(df[~df['is_similar']], columns = COLUMNS)
    elif output_option == 'all':
        df_sim = pd.DataFrame(df, columns = COLUMNS + ['is_similar'])
    else:
        raise Exception('Error: Invalid output option.')

    # write to file
    if method in {'fingerprint', 'mcs'}:
        df_sim.sort_values(by=['Similarity_Score'], ascending=False, inplace=True)
    df_sim = df_sim.reset_index(drop=True)
    print('Number of rows:', df_sim.shape[0])
    df_sim = remove_unnamed_columns(df_sim)
    df_sim.to_csv(output_file + f'_{df_sim.shape[0]}.csv')


def similarity_search_multiple_query(input_file_library, input_file_query, id_column_name_query, smiles_column_name_query, method, **kwargs):
    """
    Find SMILES in input_file_library that are similar to the SMILES in input_file_query.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param input_file_library: str, path of the input library.
    :param input_file_query: str, path of the input query SMILES.
    :param id_column_name_query: str, name of the ID column in input_file_query.
    :param smiles_column_name_query: str, name of the SMILES column in input_file_query.
    :param method: str, method for similarity search, allowed values include 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param mcs_match: str, specify 'mcs_match' to use different comparison functions for MCS, allowed values include 'exact', 'anyAtom'.
    :param substructure_method: str, method used to convert smiles_query to mol_query, allowed values include 'SMILES', 'SMARTS'.
    :param output_folder: str, path of the output file.
    :param output_file: str, name of the output file.
    :param output_option: str, options to output data, allowed values include 'satisfied', 'not_satisfied' and 'all'.
    """
    # read query SMILES
    df_query = pd.read_csv(input_file_query)
    IDs = [str(id) for id in df_query[id_column_name_query].values.tolist()]
    SMILES = df_query[smiles_column_name_query].values.tolist()
    num_SMILES = len(IDs)

    # run similarity search for each query compounds
    total_time = []
    for i in tqdm(range(num_SMILES), desc='Processing'):
        time.sleep(0.0001)
        start_time = time.time()

        smiles_column_name = 'Cleaned_SMILES'
        smiles_query, id = SMILES[i], IDs[i]
        print(f'Start: {id} {smiles_query}')
        output_file = f'{id}.csv'
        similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method, output_file=output_file, **kwargs)

        end_time = time.time()
        elapsed_time = end_time - start_time
        total_time.append(elapsed_time)
        print(f'Done, time = {elapsed_time:.2f}')

    print(f'Mean time is {np.mean(np.array(total_time)):.4f}, std is {np.std(np.array(total_time)):.4f}')


### Select analogs ###
def select_analogs(input_file_query, analogs_dir, analog_method, **kwargs):
    """
    Select analogs based on two criteria: 1. similarity scores are above the given similarity score cutoff;
    or 2. the rank n most similar analogs
    :param input_file_query: str, path of the input query SMILES.
    :param analogs_dir: str, path of the directory containing analogs.
    :param analog_method: str, method for selecting analogs, allowed values include 'cutoff' and 'rank'.
    :param similarity_cutoff: float, similarity cutoff for selecting analogs.
    :param similarity_rank: int, rank of the selected analog.
    :param deduplication: bool, whether to deduplicate the selected analogs
    """
    # output file
    output_file = os.path.splitext(os.path.abspath(input_file_query))[0]

    # get query IDs and SMILES
    df_query = pd.read_csv(input_file_query)
    IDs_query = df_query['ID'].tolist()
    SMILES_query = df_query['SMILES'].tolist()

    # get analogs
    if analog_method == 'cutoff':
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.5)
        print(f'Get analogs whose similarity scores are above {similarity_cutoff}.')
        output_file = f'{output_file}_cutoff{similarity_cutoff}'

        dfs = []
        files = os.listdir(analogs_dir)
        for i, id in enumerate(IDs_query):
            analogs_file = [file for file in files if file.startswith(f'{id}_')][0]
            df_analogs = pd.read_csv(f'{analogs_dir}/{analogs_file}')
            df_analogs['ID'] = df_analogs['ID'].apply(lambda analog_id: str(analog_id))

            df_selected_analogs = df_analogs[df_analogs['Similarity_Score'] >= similarity_cutoff]

            if df_selected_analogs.empty:
                continue
            df_selected_analogs = pd.DataFrame(df_selected_analogs, columns=['ID', 'SMILES', 'Similarity_Score'])
            df_selected_analogs.rename(columns={'ID': 'Analog_ID', 'SMILES': 'Analog_SMILES'}, inplace=True)
            df_selected_analogs['Query_ID'] = str(id)
            df_selected_analogs['Query_SMILES'] = SMILES_query[i]
            dfs.append(df_selected_analogs)

    elif analog_method == 'rank':
        similarity_rank = kwargs.get('similarity_rank', 1)
        assert similarity_rank >= 1, 'Error: Invalid rank.'
        print(f'Get the rank {similarity_rank} most similar analog.')
        output_file = f'{output_file}_rank{similarity_rank}'

        dfs = []
        files = os.listdir(analogs_dir)
        for i, id in enumerate(IDs_query):
            analogs_file = [file for file in files if file.startswith(f'{id}_')][0]
            df_analogs = pd.read_csv(f'{analogs_dir}/{analogs_file}')
            df_analogs['ID'] = df_analogs['ID'].apply(lambda analog_id: str(analog_id))

            df_selected_analogs = pd.DataFrame()
            if df_analogs.shape[0] >= similarity_rank:
                df_analogs = df_analogs.sort_values(by=['Similarity_Score'], ascending=[False], ignore_index=True)
                df_selected_analogs = df_analogs.loc[[similarity_rank - 1]]

            if df_selected_analogs.empty:
                continue
            df_selected_analogs = pd.DataFrame(df_selected_analogs, columns=['ID', 'SMILES', 'Similarity_Score'])
            df_selected_analogs.rename(columns={'ID': 'Analog_ID', 'SMILES': 'Analog_SMILES'}, inplace=True)
            df_selected_analogs['Query_ID'] = str(id)
            df_selected_analogs['Query_SMILES'] = SMILES_query[i]
            dfs.append(df_selected_analogs)

    else:
        raise Exception('Error: Invalid method for selecting analogs.')

    # concat results
    df_concat = pd.concat(dfs, ignore_index=True, sort=False)

    # remove duplicates
    deduplication = kwargs.get('deduplication', False)
    if deduplication:
        df_concat = df_concat.sort_values(by=['Similarity_Score'], ascending=[False], ignore_index=True)
        df_concat = df_concat.drop_duplicates('Analog_ID', keep='first', ignore_index=True)

    # write output file
    df_concat = df_concat.reset_index(drop=True)
    print('Number of selected analogs:', df_concat.shape[0])
    df_concat = remove_unnamed_columns(df_concat)
    df_concat.to_csv(f'{output_file}_{df_concat.shape[0]}.csv')


### Plot similarity score distribution ###
def plot_distribution(input_file):
    """
    Plot similarity score distribution
    :param input_file: str, path of the input file
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file))[0]
    df = pd.read_csv(input_file)

    # plot distribution
    values = df['Similarity_Score'].tolist()
    output_file = f'{output_file}_similarity_distribution.pdf'

    plt.figure(1)
    plt.hist(values, 10, range=(0.0, 1.0))
    plt.xlabel('Similarity Score', fontproperties=font)
    plt.ylabel('Counts', fontproperties=font)
    plt.xticks(fontproperties=font)
    plt.yticks(fontproperties=font)

    plt.savefig(output_file, dpi=300)
    plt.show()
    plt.close()



if __name__ == '__main__':

    single_query = False

    ### Similarity search for a single query SMILES ###
    if single_query:
        input_file_library = 'similarity_search/tests/library.csv'
        # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
        smiles_column_name = 'Cleaned_SMILES'
        smiles_query = 'N#CC1(Oc2ccc(-c3ccccc3)cc2)CC2(CCN(C(=O)Nc3ccccc3)C2)C1'
        # smiles_query = 'C1=CC=CC=C1'

        output_folder = 'similarity_search/tests/similarity_search_results'
        output_file = 'output.csv'
        output_option = 'satisfied'

        ### Method: 'fingerprint' ###
        # similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method='fingerprint', similarity_cutoff=0.3, output_folder=output_folder, output_file=output_file, output_option=output_option)
        ### Method: 'mcs' ###
        # similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method='mcs', mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_file=output_file, output_option=output_option)
        ### Method: 'substructure' ###
        # similarity_search_single_query(input_file_library, smiles_column_name, smiles_query, method='substructure', substructure_method='SMARTS', output_folder=output_folder, output_file=output_file, output_option=output_option)

    ### Similarity search for multiple SMILES ###
    else:
        input_file_library = 'similarity_search/tests/library.csv'
        # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
        input_file_query = 'similarity_search/tests/query_cmps.csv'
        id_column_name_query = 'ID'
        smiles_column_name_query = 'Cleaned_SMILES'

        output_folder = 'similarity_search/tests/similarity_search_results'
        output_option = 'satisfied'

        ### Method: 'fingerprint' ###
        # similarity_search_multiple_query(input_file_library, input_file_query, id_column_name_query, smiles_column_name_query,
        #                                method='fingerprint', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
        ### Method: 'mcs' ###
        # similarity_search_multiple_query(input_file_library, input_file_query, id_column_name_query, smiles_column_name_query,
        #                                method='mcs', mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
        ### Method: 'substructure' ###
        # similarity_search_multiple_query(input_file_library, input_file_query, id_column_name_query, smiles_column_name_query,
        #                                method='substructure', substructure_method='SMARTS',output_folder=output_folder, output_option=output_option)


    ### Select analogs ###
    input_file_query = 'similarity_search/tests/query_cmps.csv'
    analogs_dir = 'similarity_search/tests/similarity_search_results'
    ### Select analogs above similarity cutoff
    # analog_method = 'cutoff'
    # similarity_cutoff = 0.4
    # select_analogs(input_file_query, analogs_dir, analog_method, similarity_cutoff=similarity_cutoff, deduplication=False)
    ### Select the rank-n most similar analog
    # analog_method = 'rank'
    # similarity_rank = 1
    # select_analogs(input_file_query, analogs_dir, analog_method, similarity_rank=similarity_rank, deduplication=False)


    ### Plot similarity score distribution ###
    # input_file = 'similarity_search/tests/query_cmps_rank1_5.csv'
    # plot_distribution(input_file)



    ### Tests ###
    # smiles = 'C1=CC=CC=C1'
    # smiles_query = 'C1=CC=CC=C1'   # True

    # smiles = 'Clc1ccc(CC)cc1'   # True
    # smiles = 'Clc1ccc(CC)c(C)c1'   # True
    # smiles = 'c1ccccc1'   # False
    # smiles_query = '[*]c1ccc([*])cc1'

    # smiles = 'Clc1c([H])c([H])c(CC)c([H])c1[H]'   # True
    # smiles = 'Clc1c([H])c([H])c(CC)c(C)c1[H]'   # True
    # smiles = '[H]c1c([H])c([H])c([H])c([H])c1[H]'   # False
    # smiles_query = '[*]c1c([H])c([H])c([*])c([H])c1[H]'

    # smiles = 'NCC'   # True
    # smiles_query = '[*]N'

    # smiles = 'C1CC[NH2+]CC1'   # True
    # smiles = '[NH3+]CC'   # True
    # smiles_query = '[*][N+]'

    # smiles = 'O=C(C1=CC(N2CCNCC2)=CC=C1C)NC3(CC3)c4c5ccccc5ccc4'
    # smiles = 'O=C(C1=CC(N2CCNCC2)=CC=C1C)NCCc3c4ccccc4ccc3'

    # smiles_query = '[*]C1=CC=CC2=C1C=CC=C2'   # False if using kekulized form
    # smiles_query = '*c1cccc2ccccc12'
    # smiles_query = '[*]NC([*])=O'
    # smiles_query = '*NC(*)=O'
    # smiles_query = 'O=C([*]N)N[*]c1c2ccccc2ccc1'

    # smiles_query = '[*]C(N[*]c(ccc1)c2c1cccc2)=O'   # True
    # smiles_query = 'O=C(c1cccc(N)c1)N[*]c(ccc2)c3c2cccc3'   # True
    # smiles_query = 'O=C(c1cccc([*]N)c1)N[*]c(ccc2)c3c2cccc3'
    # smiles_query = 'c1(cccc2ccccc12)*NC(=O)*'

    # print(is_similar_by_substructure(smiles, smiles_query, substructure_method='SMARTS'))







