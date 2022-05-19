#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 10:50:25 2022

@author: guohan
"""

import os
import warnings
import sys

path_list = sys.path
module_path = '/Users/guohan/Documents/Code/Tool/utils'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

import pandas as pd
import numpy as np

from molecular_description import get_fingerprint


if __name__ == '__main__':

    input_file = 'TRPML1_activity_inhousePatentPubchem_class2_subset4_96772.csv'
    df = pd.read_csv(input_file, index_col = 0)
    
    IDs = df['ID'].values.tolist()
    SMILES = df['Cleaned_SMILES'].values.tolist()
    
    for i in range(len(IDs)):
        print(IDs[i])
        get_fingerprint(SMILES[i], fp_method = 'ecfp4')
    


