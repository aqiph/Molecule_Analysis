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
def is_similar_by_fingerprint(smiles, smiles_ref, similarity_cutoff):
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
    return similarity_score >= similarity_cutoff


def is_similar_by_mcs(smiles, smiles_ref, similarity_cutoff):
    """
    Determine whether smiles is similar to smiles_ref based on MCS
    :param smiles: str, SMILES
    :param smiles_ref: str, reference SMILES
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles is considered to be similar to smiles_ref
    :return: bool, whether smiles is considered to be similar to smiles_ref
    """
    similarity_score = cal_MCS(smiles, smiles_ref)
    return similarity_score >= similarity_cutoff


def is_similar_by_substructure(smiles, smiles_ref):
    """
    Determine whether smiles_ref is present in smiles
    :param smiles: str, SMILES
    :param smiles_ref: str, reference SMILES
    :return: bool, whether smiles is considered to be similar to smiles_ref
    """
    # check reference SMILES
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
    mol_ref = Chem.MolFromSmiles(smiles_ref)
    return mol.HasSubstructMatch(mol_ref)


def similarity_search_single_ref(df, smiles_column_name, smiles_ref, method, similarity_cutoff=0.8):
    """
    Determine whether SMILES in df is similar to the reference SMILES.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param df: pd.DataFrame object, library whose smiles_column_name column contains compound SMILES to be searched from.
    :param smiles_column_name: str, name of the SMILES column
    :param smiles_ref: str, reference SMILES
    :param method: str, method for similarity search, allowed values contain 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
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
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_fingerprint(smiles, smiles_ref, similarity_cutoff))
    elif method == 'mcs':
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_mcs(smiles, smiles_ref, similarity_cutoff))
    elif method == 'substructure':
        df['is_similar'] = df[smiles_column_name].apply(lambda smiles: is_similar_by_substructure(smiles, smiles_ref))
    else:
        raise Exception('Error: Invalid similarity search method.')

    return df


def main(input_file, smiles_column_name, smiles_ref, method, similarity_cutoff=0.8, output_file = None, output_option='selected'):
    """
    Determine whether SMILES in input file is similar to the reference SMILES.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param input_file: str, path of the input library
    :param smiles_column_name: str, name of the SMILES column
    :param smiles_ref: str, reference SMILES
    :param method: str, method for similarity search, allowed values contain 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param output_file: str, name of the output file
    :param output_option: str, options to output data, allowed values include 'selected', 'not_selected' and 'all'
    """
    # read df
    df = pd.read_csv(input_file)
    COLUMNS = df.columns.tolist()
    assert smiles_column_name in COLUMNS, 'Error: Invalid input SMILES column name.'

    # output file
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:
        output_file = basename
    output_file = os.path.join(folder, os.path.splitext(output_file)[0])

    # similarity search
    df = similarity_search_single_ref(df, smiles_column_name, smiles_ref, method, similarity_cutoff)
    if output_option == 'selected':
        df_sim = pd.DataFrame(df[df['is_similar']], columns = COLUMNS)
    elif output_option == 'not_selected':
        df_sim = pd.DataFrame(df[~df['is_similar']], columns = COLUMNS)
    elif output_option == 'all':
        df_sim = pd.DataFrame(df, columns = COLUMNS + ['is_similar'])
    else:
        raise Exception('Error: Invalid output option.')

    # compute mcs score
    df_sim['MCS_Score'] = df_sim[smiles_column_name].apply(lambda smiles: cal_MCS(smiles, smiles_ref))

    # write to file
    df_sim.sort_values(by=['MCS_Score'], ascending=False, inplace=True)
    df_sim = df_sim.reset_index(drop=True)
    print('Number of rows:', df_sim.shape[0])
    df_sim = remove_unnamed_columns(df_sim)
    df_sim.to_csv(output_file + f'_{df_sim.shape[0]}.csv')



if __name__ == '__main__':

    # input_file = 'similarity_search/tests/library.csv'
    input_file = 'similarity_search/tests/HTS_forGNN_446663.csv'
    # input_file = 'similarity_search/tests/HTS_Expresspick-Pfizer_forGNN_4208.csv'
    smiles_column_name = 'Cleaned_SMILES'
    smiles_ref = 'NC(=O)SCC(=O)Nc1ccc(S(N)(=O)=O)cc1'
    # smiles_ref = 'O=C(NCc1ccc(F)cc1)c1cccc(Nc2ccc(-c3ccccc3)nn2)c1'
    method = 'mcs'
    output_file = 'output.csv'
    main(input_file, smiles_column_name, smiles_ref, method, similarity_cutoff=0.6, output_file=output_file, output_option='selected')

