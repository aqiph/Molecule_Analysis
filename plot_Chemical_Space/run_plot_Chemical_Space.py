#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 07 16:36:00 2024

@author: guohan
"""

import sys
path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Molecule_Analysis'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

from plot_Chemical_Space import plot_chemical_space_tSNE

if __name__=='__main__':
    input_file_library = 'tests/test_plot_Chemical_Space_train.csv'
    input_file_query = 'tests/test_plot_Chemical_Space_test.csv'
    plot_chemical_space_tSNE(input_file_library,
                             input_file_query=input_file_query,
                             smiles_column_name_library='SMILES',
                             smiles_column_name_query='SMILES',
                             label_column_name_library=None,
                             label_column_name_query=None,
                             plot_order=None)