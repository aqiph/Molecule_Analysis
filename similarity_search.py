#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:10:00 2023

@author: guohan

1. Similarity search
2. Get most similar analogs
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
def is_similar_by_fingerprint(smiles, smiles_ref, similarity_cutoff=0.7):
    """
    Determine whether smiles is similar to smiles_ref based on ecfp4 fingerprint
    :param smiles: str, SMILES
    :param smiles_ref: str, reference SMILES
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles is considered to be similar to smiles_ref
    :return: bool, whether smiles is considered to be similar to smiles_ref
    """
    fp_1 = get_fingerprint(smiles, fp_method = 'ecfp4')
    fp_2 = get_fingerprint(smiles_ref, fp_method = 'ecfp4')
    dist = cal_fingerprint_distance(fp_1, fp_2, fp_method='ecfp4')
    similarity_score = 1.0 - dist
    return similarity_score >= similarity_cutoff, similarity_score


def is_similar_by_mcs(smiles, smiles_ref, match='exact', similarity_cutoff=0.7):
    """
    Determine whether smiles is similar to smiles_ref based on MCS
    :param smiles: str, SMILES
    :param smiles_ref: str, reference SMILES
    :param match: str, specify 'match' to use different comparison functions, allowed values include 'exact', 'anyAtom'
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles is considered to be similar to smiles_ref
    :return: bool, whether smiles is considered to be similar to smiles_ref
    """
    similarity_score = cal_MCS(smiles, smiles_ref, match)
    return similarity_score >= similarity_cutoff, similarity_score


def is_similar_by_substructure(smiles, smiles_ref, substructure_method='SMILES'):
    """
    Determine whether smiles_ref is present in smiles
    :param smiles: str, SMILES
    :param smiles_ref: str, reference SMILES
    :param substructure_method: str, method used to convert smiles_ref to mol_ref, allowed values include 'SMILES', 'SMARTS'
    :return: bool, whether smiles is considered to be similar to smiles_ref
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

    # determine whether smiles_ref is present in smiles
    if substructure_method == 'SMILES':
        mol_ref = Chem.MolFromSmiles(smiles_ref)
    elif substructure_method == 'SMARTS':
        smiles_ref = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_ref))   # Kekulize molecule before get mol as a pattern
        # print(smiles_ref)
        mol_ref = Chem.MolFromSmarts(smiles_ref)
    else:
        raise Exception('Error: Invalid substructure method.')

    hasSubstructure = mol.HasSubstructMatch(mol_ref)
    return hasSubstructure, int(hasSubstructure)


def similarity_search_single_ref(input_file, smiles_column_name, smiles_ref, method, **kwargs):
    """
    Find SMILES in input_file that are similar to the reference SMILES.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param input_file: str, path of the input library
    :param smiles_column_name: str, name of the SMILES column
    :param smiles_ref: str, reference SMILES
    :param method: str, method for similarity search, allowed values include 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param mcs_match: str, specify 'mcs_match' to use different comparison functions for MCS, allowed values include 'exact', 'anyAtom'.
    :param substructure_method: str, method used to convert smiles_ref to mol_ref, allowed values include 'SMILES', 'SMARTS'
    :param output_folder: str, path of the output file
    :param output_file: str, name of the output file
    :param output_option: str, options to output data, allowed values include 'satisfied', 'not_satisfied' and 'all'
    """
    # files
    df = pd.read_csv(input_file)
    COLUMNS = df.columns.tolist() + ['Similarity_Score']
    assert smiles_column_name in COLUMNS, 'Error: Invalid input SMILES column name.'

    folder = kwargs.get('output_folder', os.getcwd())
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f'Fold {folder} created.')
    output_file = kwargs.get('output_file', input_file)
    output_file = os.path.splitext(os.path.split(output_file)[1])[0]
    output_file = os.path.join(folder, output_file)

    # check reference SMILES
    if not smiles_ref:
        print(f"Error: Invalid SMILES {smiles_ref}")
        return None
    try:
        mol = Chem.MolFromSmiles(smiles_ref)
    except Exception:
        print(f"Error: Invalid SMILES {smiles_ref}")
        return None
    if not mol:
        print(f"Error: Invalid SMILES {smiles_ref}")
        return None

    # determine whether SMILES in df is similar to the reference SMILES
    method = method.lower()
    if method == 'fingerprint':
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_fingerprint(smiles, smiles_ref, similarity_cutoff))
    elif method == 'mcs':
        match = kwargs.get('mcs_match', 'exact')
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_mcs(smiles, smiles_ref, match, similarity_cutoff))
    elif method == 'substructure':
        substructure_method = kwargs.get('substructure_method', 'SMILES')
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_substructure(smiles, smiles_ref, substructure_method))
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


def similarity_search_multiple_ref(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref, method, **kwargs):
    """
    Find SMILES in input_file_library that are similar to the SMILES in input_file_ref.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param input_file_library: str, path of the input library.
    :param input_file_ref: str, path of the input reference SMILES.
    :param id_column_name_ref: str, name of the ID column in input_file_ref
    :param smiles_column_name_ref: str, name of the SMILES column in input_file_ref
    :param method: str, method for similarity search, allowed values include 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param mcs_match: str, specify 'mcs_match' to use different comparison functions for MCS, allowed values include 'exact', 'anyAtom'.
    :param substructure_method: str, method used to convert smiles_ref to mol_ref, allowed values include 'SMILES', 'SMARTS'
    :param output_folder: str, path of the output file
    :param output_file: str, name of the output file
    :param output_option: str, options to output data, allowed values include 'satisfied', 'not_satisfied' and 'all'
    """
    # read reference SMILES
    df_ref = pd.read_csv(input_file_ref)
    IDs = df_ref[id_column_name_ref].values.tolist()
    SMILES = df_ref[smiles_column_name_ref].values.tolist()
    num_SMILES = len(IDs)

    # run similarity search for each reference cmps
    total_time = []
    for i in tqdm(range(num_SMILES), desc='Processing'):
        time.sleep(0.0001)
        start_time = time.time()

        smiles_column_name = 'Cleaned_SMILES'
        smiles_ref, id = SMILES[i], IDs[i]
        print(f'Start: {id} {smiles_ref}')
        output_file = f'{id}.csv'
        similarity_search_single_ref(input_file_library, smiles_column_name, smiles_ref, method, output_file=output_file, **kwargs)

        end_time = time.time()
        elapsed_time = end_time - start_time
        total_time.append(elapsed_time)
        print(f'Done, time = {elapsed_time:.2f}')

    print(f'Mean time is {np.mean(np.array(total_time)):.4f}, std is {np.std(np.array(total_time)):.4f}')


### Get most similar analogs ###
def get_most_similar_analog(input_file, analogs_dir, n=1):
    """
    Get the rank-n most similar analogs
    :param input_file: str, path of the file containing reference compounds
    :param analogs_dir: str, path of the folder containing analogs
    :param n: int, specify the rank of the similar analog
    """
    assert n >= 1, 'Error: Invalid rank.'

    # output file
    output_file = os.path.splitext(os.path.abspath(input_file))[0]

    df = pd.read_csv(input_file)
    df_bestScore = pd.DataFrame(df, columns = ['ID', 'SMILES', 'Cleaned_SMILES'])

    df_bestScore['Analog'] = df_bestScore['ID'].apply(lambda cmp: get_analog(cmp, analogs_dir, n-1))
    df_bestScore[['Analog_ID', 'Analog_SMILES', 'Analog_Cleaned_SMILES', 'Similarity_Score']] = \
        pd.DataFrame(df_bestScore['Analog'].tolist())
    COLUMNS = ['ID', 'SMILES', 'Cleaned_SMILES', 'Analog_ID', 'Analog_SMILES', 'Analog_Cleaned_SMILES', 'Similarity_Score']
    df_bestScore = pd.DataFrame(df_bestScore, columns = COLUMNS)

    df_bestScore = df_bestScore.reset_index(drop=True)
    print('Number of rows:', df_bestScore.shape[0])
    df_bestScore = remove_unnamed_columns(df_bestScore)
    df_bestScore.to_csv(f'{output_file}_Top{n}_{df_bestScore.shape[0]}.csv')


def get_analog(cmp, analogs_dir, n):
    files = os.listdir(analogs_dir)
    filtered_files = [file for file in files if file.startswith(f'{cmp}_')]
    df_cmp = pd.read_csv(f'{analogs_dir}/{filtered_files[0]}')

    if df_cmp.shape[0] <= n:
        return ['', '', '', 0.0]

    return df_cmp.iloc[n, [False, True, True, True, False, False, False, False, False, True]].tolist()


def plot_distribution(input_file):
    """
    Plot similarity score distribution
    """
    # files
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

    single_ref = False

    ### Similarity search for a single reference SMILES ###
    if single_ref:
        input_file = 'similarity_search/tests/library.csv'
        # input_file = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
        smiles_column_name = 'Cleaned_SMILES'
        smiles_ref = 'N#CC1(Oc2ccc(-c3ccccc3)cc2)CC2(CCN(C(=O)Nc3ccccc3)C2)C1'
        # smiles_ref = 'C1=CC=CC=C1'

        output_folder = 'similarity_search/tests/similarity_search_results'
        output_file = 'output.csv'
        output_option = 'satisfied'

        ### Method: 'fingerprint' ###
        # similarity_search_single_ref(input_file, smiles_column_name, smiles_ref, method='fingerprint', similarity_cutoff=0.3, output_folder=output_folder, output_file=output_file, output_option=output_option)
        ### Method: 'mcs' ###
        # similarity_search_single_ref(input_file, smiles_column_name, smiles_ref, method='mcs', mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_file=output_file, output_option=output_option)
        ### Method: 'substructure' ###
        # similarity_search_single_ref(input_file, smiles_column_name, smiles_ref, method='substructure', substructure_method='SMARTS', output_folder=output_folder, output_file=output_file, output_option=output_option)

    ### Similarity search for multiple SMILES ###
    else:
        input_file_library = 'similarity_search/tests/library.csv'
        # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
        input_file_ref = 'similarity_search/tests/reference_cmps.csv'
        id_column_name_ref = 'ID'
        smiles_column_name_ref = 'Cleaned_SMILES'

        output_folder = 'similarity_search/tests/similarity_search_results'
        output_option = 'satisfied'

        ### Method: 'fingerprint' ###
        # similarity_search_multiple_ref(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref,
        #                                method='fingerprint', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
        ### Method: 'mcs' ###
        # similarity_search_multiple_ref(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref,
        #                                method='mcs', mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
        ### Method: 'substructure' ###
        # similarity_search_multiple_ref(input_file_library, input_file_ref, id_column_name_ref, smiles_column_name_ref,
        #                                method='substructure', substructure_method='SMARTS',output_folder=output_folder, output_option=output_option)


    ### Get the rank-n most similar analogs ###
    # input_file = 'similarity_search/tests/reference_cmps.csv'
    # analogs_dir = 'similarity_search/tests/similarity_search_results'
    # get_most_similar_analog(input_file, analogs_dir, n=1)


    ### Plot similarity score distribution
    # input_file = 'similarity_search/tests/reference_cmps_Top1_12.csv'
    # plot_distribution(input_file)



    ### Tests ###
    # smiles = 'C1=CC=CC=C1'
    # smiles_ref = 'C1=CC=CC=C1'   # True

    # smiles = 'Clc1ccc(CC)cc1'   # True
    # smiles = 'Clc1ccc(CC)c(C)c1'   # True
    # smiles = 'c1ccccc1'   # False
    # smiles_ref = '[*]c1ccc([*])cc1'

    # smiles = 'Clc1c([H])c([H])c(CC)c([H])c1[H]'   # True
    # smiles = 'Clc1c([H])c([H])c(CC)c(C)c1[H]'   # True
    # smiles = '[H]c1c([H])c([H])c([H])c([H])c1[H]'   # False
    # smiles_ref = '[*]c1c([H])c([H])c([*])c([H])c1[H]'

    # smiles = 'NCC'   # True
    # smiles_ref = '[*]N'

    # smiles = 'C1CC[NH2+]CC1'   # True
    # smiles = '[NH3+]CC'   # True
    # smiles_ref = '[*][N+]'

    # smiles = 'O=C(C1=CC(N2CCNCC2)=CC=C1C)NC3(CC3)c4c5ccccc5ccc4'
    # smiles = 'O=C(C1=CC(N2CCNCC2)=CC=C1C)NCCc3c4ccccc4ccc3'

    # smiles_ref = '[*]C1=CC=CC2=C1C=CC=C2'   # False if not kekulized
    # smiles_ref = '*c1cccc2ccccc12'
    # smiles_ref = '[*]NC([*])=O'
    # smiles_ref = '*NC(*)=O'
    # smiles_ref = 'O=C([*]N)N[*]c1c2ccccc2ccc1'

    # smiles_ref = '[*]C(N[*]c(ccc1)c2c1cccc2)=O'   # True
    # smiles_ref = 'O=C(c1cccc(N)c1)N[*]c(ccc2)c3c2cccc3'   # True
    # smiles_ref = 'O=C(c1cccc([*]N)c1)N[*]c(ccc2)c3c2cccc3'
    # smiles_ref = 'c1(cccc2ccccc12)*NC(=O)*'



    # print(is_similar_by_substructure(smiles, smiles_ref, substructure_method='SMARTS'))







