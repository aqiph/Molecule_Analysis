#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:10:00 2023

@author: guohan

1. Similarity search

"""


import os
import pandas as pd
from rdkit import Chem

from utils.molecular_description import get_fingerprint, cal_fingerprint_distance, cal_MCS
from utils.tools import remove_unnamed_columns



### Helper function ###
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


def is_similar_by_mcs(smiles, smiles_ref, similarity_cutoff=0.7):
    """
    Determine whether smiles is similar to smiles_ref based on MCS
    :param smiles: str, SMILES
    :param smiles_ref: str, reference SMILES
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles is considered to be similar to smiles_ref
    :return: bool, whether smiles is considered to be similar to smiles_ref
    """
    similarity_score = cal_MCS(smiles, smiles_ref)
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
        print(smiles_ref)
        mol_ref = Chem.MolFromSmarts(smiles_ref)
    else:
        raise Exception('Error: Invalid substructure method.')

    hasSubstructure = mol.HasSubstructMatch(mol_ref)
    return hasSubstructure, int(hasSubstructure)


def similarity_search_single_ref(df, smiles_column_name, smiles_ref, method, **kwargs):
    """
    Determine whether SMILES in df is similar to the reference SMILES.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param df: pd.DataFrame object, library whose smiles_column_name column contains compound SMILES to be searched from.
    :param smiles_column_name: str, name of the SMILES column
    :param smiles_ref: str, reference SMILES
    :param method: str, method for similarity search, allowed values include 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param substructure_method: str, method used to convert smiles_ref to mol_ref, allowed values include 'SMILES', 'SMARTS'
    :return: df, pd.DataFrame object, library with one additional column specifying whether compounds are similar to the reference compound.
    """
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

    # Determine whether SMILES in df is similar to the reference SMILES
    method = method.lower()
    if method == 'fingerprint':
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_fingerprint(smiles, smiles_ref, similarity_cutoff))
    elif method == 'mcs':
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_mcs(smiles, smiles_ref, similarity_cutoff))
    elif method == 'substructure':
        substructure_method = kwargs.get('substructure_method', 'SMILES')
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_substructure(smiles, smiles_ref, substructure_method))
    else:
        raise Exception('Error: Invalid similarity search method.')

    df[['is_similar', 'Similarity_Score']] = pd.DataFrame(df['is_similar'].values.tolist())

    return df


def main(input_file, smiles_column_name, smiles_ref, method, output_file = None, **kwargs):
    """
    Determine whether SMILES in input file is similar to the reference SMILES.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param input_file: str, path of the input library
    :param smiles_column_name: str, name of the SMILES column
    :param smiles_ref: str, reference SMILES
    :param method: str, method for similarity search, allowed values include 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param output_file: str, name of the output file
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param substructure_method: str, method used to convert smiles_ref to mol_ref, allowed values include 'SMILES', 'SMARTS'
    :param output_option: str, options to output data, allowed values include 'selected', 'not_selected' and 'all'
    """
    # read df
    df = pd.read_csv(input_file)
    COLUMNS = df.columns.tolist() + ['Similarity_Score']
    assert smiles_column_name in COLUMNS, 'Error: Invalid input SMILES column name.'

    # output file
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:
        output_file = basename
    output_file = os.path.join(folder, os.path.splitext(output_file)[0])
    output_option = kwargs.get('output_option', 'selected')

    # similarity search
    df = similarity_search_single_ref(df, smiles_column_name, smiles_ref, method, **kwargs)
    if output_option == 'selected':
        df_sim = pd.DataFrame(df[df['is_similar']], columns = COLUMNS)
    elif output_option == 'not_selected':
        df_sim = pd.DataFrame(df[~df['is_similar']], columns = COLUMNS)
    elif output_option == 'all':
        df_sim = pd.DataFrame(df, columns = COLUMNS + ['is_similar'])
    else:
        raise Exception('Error: Invalid output option.')

    # write to file
    df_sim.sort_values(by=['Similarity_Score'], ascending=False, inplace=True)
    df_sim = df_sim.reset_index(drop=True)
    print('Number of rows:', df_sim.shape[0])
    df_sim = remove_unnamed_columns(df_sim)
    df_sim.to_csv(output_file + f'_{df_sim.shape[0]}.csv')



if __name__ == '__main__':

    # input_file = 'similarity_search/tests/library.csv'
    # input_file = 'similarity_search/tests/HTS_forGNN_446663.csv'
    # smiles_column_name = 'Cleaned_SMILES'
    # smiles_ref = 'c1ccccc1'
    # output_file = 'output.csv'
    ### Method: 'fingerprint' ###
    # main(input_file, smiles_column_name, smiles_ref, method='fingerprint', output_file=output_file, similarity_cutoff=0.3, output_option='all')
    ### Method: 'mcs' ###
    # main(input_file, smiles_column_name, smiles_ref, method='mcs', output_file=output_file, similarity_cutoff=0.3, output_option='all')
    ### Method: 'substructure' ###
    # main(input_file, smiles_column_name, smiles_ref, method='substructure', output_file=output_file, substructure_method='SMARTS', output_option='all')


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

    smiles = 'O=C(C1=CC(N2CCNCC2)=CC=C1C)NC3(CC3)c4c5ccccc5ccc4'
    # smiles = 'O=C(C1=CC(N2CCNCC2)=CC=C1C)NCCc3c4ccccc4ccc3'

    # smiles_ref = '[*]C1=CC=CC2=C1C=CC=C2'   # False if not kekulized
    # smiles_ref = '*c1cccc2ccccc12'
    # smiles_ref = '[*]NC([*])=O'
    # smiles_ref = '*NC(*)=O'
    # smiles_ref = 'O=C([*]N)N[*]c1c2ccccc2ccc1'

    # smiles_ref = '[*]C(N[*]c(ccc1)c2c1cccc2)=O'   # True
    # smiles_ref = 'O=C(c1cccc(N)c1)N[*]c(ccc2)c3c2cccc3'   # True
    # smiles_ref = 'O=C(c1cccc([*]N)c1)N[*]c(ccc2)c3c2cccc3'
    smiles_ref = 'c1(cccc2ccccc12)*NC(=O)*'



    print(is_similar_by_substructure(smiles, smiles_ref, substructure_method='SMARTS'))







