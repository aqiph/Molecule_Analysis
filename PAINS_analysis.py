#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 06 16:58:00 2024

@author: guohan

"""


import os
import pandas as pd
from rdkit import Chem
from utils.tools import remove_unnamed_columns


class PAINS_Filter(object):
    def __init__(self):
        self._get_pains_substructures()


    def _get_pains_substructures(self):
        """
        Private function for reading pre-defined PAINS substructures from SMARTS
        """
        pains_file = '/Users/guohan/Documents/Codes/Molecule_Analysis/utils/pains.txt'
        with open(pains_file) as fh_pains_substructure_SMARTS:
            pains_substructure_SMARTS = fh_pains_substructure_SMARTS.readlines()
            self.pains_substructure_list = [Chem.MolFromSmarts(smarts.strip()) for smarts in pains_substructure_SMARTS]
            print(f'Get {len(self.pains_substructure_list)} PAINS substructures.')


    def _pains_filter(self, mol):
        """
        Private function for checking whether a compound in mol object contains PAINS substructures.
        :param mol: mol object, query compound
        """
        for pains_substructure in self.pains_substructure_list:
            if mol.HasSubstructMatch(pains_substructure):
                return True
        return False


    def is_PAINS(self, smiles):
        """
        Check whether a SMILES is PAINS.
        Refer to Baell and Holloway, J. Med. Chem. 2010, 53(7), 2719-2740,
        and Stephen J. et al., J. Chem. Inf. Model., 2017, 57, 417−427.
        :param smiles: str, query compound.
        """
        # check SMILES
        if not smiles:
            print(f"Error: Invalid SMILES {smiles}")
            return True
        try:
            mol = Chem.MolFromSmiles(smiles)
        except Exception:
            print(f"Error: Invalid SMILES {smiles}")
            return True
        if not mol:
            print(f"Error: Invalid SMILES {smiles}")
            return True

        return self._pains_filter(mol)


    def filter_PAINS(self, input_file, smiles_column_name, output_option='notPAINS'):
        """
        Recognize PAINS in input_file based on substructures and output cleaned SMILES, PAINS or PAINS labels based on
        output option.
        Refer to Baell and Holloway, J. Med. Chem. 2010, 53(7), 2719-2740,
        and Stephen J. et al., J. Chem. Inf. Model., 2017, 57, 417−427.
        :param input_file: str, path of the input file.
        :param smiles_column_name: str, name of the SMILES column in input_file.
        :param output_option: str, options to output data, allowed values include 'notPAINS', 'PAINS' and 'all'.
        """
        # files
        output_file = os.path.splitext(os.path.abspath(input_file))[0] + f'_PAINS-Filter-{output_option}'
        df = pd.read_csv(input_file)
        COLUMNS = df.columns.tolist()

        # filter PAINS
        df['PAINS_Label'] = df[smiles_column_name].apply(lambda smiles: self.is_PAINS(smiles))
        if output_option == 'notPAINS':
            df = pd.DataFrame(df[~df['PAINS_Label']], columns=COLUMNS)
        elif output_option == 'PAINS':
            df = pd.DataFrame(df[df['PAINS_Label']], columns=COLUMNS)
        elif output_option == 'all':
            df['PAINS_Label'] = df['PAINS_Label'].apply(lambda b: 'PAINS' if b else 'notPAINS')
            df = pd.DataFrame(df, columns = COLUMNS + ['PAINS_Label'])
        else:
            raise Exception('Error: Invalid output option.')

        # write to file
        df = df.reset_index(drop=True)
        print('Number of compounds:', df.shape[0])
        df = remove_unnamed_columns(df)
        df.to_csv(f'{output_file}_{df.shape[0]}.csv')



if __name__ == '__main__':

    PAINS_filter = PAINS_Filter()
    input_file = 'PAINS_analysis/tests/test_PAINS_analysis.csv'
    smiles_column_name = 'SMILES'
    PAINS_filter.filter_PAINS(input_file, smiles_column_name, output_option='all')
