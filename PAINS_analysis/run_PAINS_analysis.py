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

from PAINS_analysis import PAINS_Filter

if __name__=='__main__':
    PAINS_filter = PAINS_Filter()
    input_file = 'tests/test_PAINS_analysis.csv'
    # input_file = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGeneralUse/HTS_forGeneralUse_446664.csv'
    smiles_column_name = 'Cleaned_SMILES'
    PAINS_filter.filter_PAINS(input_file, smiles_column_name, output_option='PAINS')