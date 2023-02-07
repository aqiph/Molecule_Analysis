#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 13:22:02 2021

@author: guohan
"""
import os, sys
import shutil
import numpy as np
import pandas as pd

path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Molecule_Analysis'
if module_path not in sys.path:
    sys.path.append(module_path)
    sys.path.append(os.path.join(module_path, 'utils'))
    print('Add module path')

from SMILES_to_3D import Molecule, process_multiple_smiles
from util import cleanup_smiles


############################################
### Randomly generate SMILES from a file ###
############################################


# randomly generate a list of smiles with their IDs
def rand_smiles_list(smiles_filename, num_smiles, smiles_col_name, id_col_name):
    """
    randomly choose num_smiles of smiles from input_file_name
    :param smiles_filename:
    :param num_smiles:
    :param smiles_col_name: str, the column name of the smiles
    :param id_col_name: str, the column name of the id
    :return:
    """
    # file name
    folder, basename = os.path.split(os.path.abspath(smiles_filename))
    fmt = os.path.splitext(basename)[1]
    
    # get DataFrame from input file
    if fmt in {'.xls', '.xlsx'}:
        df = pd.read_excel(smiles_filename)
    elif fmt == '.csv':
        df = df = pd.read_csv(smiles_filename)
    else:
        print('smiles input file is not readed')
        return
    
    sample_num = df.shape[0]
    
    # get smiles and IDs
    if smiles_col_name not in df:
        raise ValueError('SMILES column name "%s" not found in csv file' % smiles_col_name)
    else:
        SMILES_list = df[smiles_col_name].to_numpy()
    
    if id_col_name not in df:
        raise ValueError('ID column name "%s" not found in csv file' % id_col_name)
    else:
        IDs = df[id_col_name].to_numpy()
    
    # shuffle data set
    shuffle_idxs = np.arange(sample_num)
    np.random.shuffle(shuffle_idxs)
    SMILES_list = SMILES_list[shuffle_idxs]
    IDs = IDs[shuffle_idxs]
    
    # get num_smiles of smiles
    SMILES_list = SMILES_list[: num_smiles]
    IDs = IDs[: num_smiles]
    
    with open(os.path.join(folder, 'smiles_list.txt'), 'w') as smiles_list:
        smiles_list.write('ID   original IDs  SMILES \n')
        for i in range(num_smiles):
            outstr = '{}   {}   '.format(i, IDs[i]) + SMILES_list[i] + '\n'
            smiles_list.write(outstr)
    
    return SMILES_list, IDs



# convert a smiles to xyz input
def smiles_to_xyz(name, smiles, cleanup, opt, plot2D):
    """
    convert a smiles to a .xyz file
    :param name: str
    :param smiles: str, smiles string
    :param cleanup: bool, whether to clean up smiles using chembl_structure_pipeline
    :return: None
    """
    assert type(smiles) == str and len(smiles) > 0
    
    # cleanup smiles
    if cleanup:
        smiles = cleanup_smiles(smiles, cleanup_chirality = False)
        
    # generate Molecuule object and 3D structure    
    molecule = Molecule(name, smiles)
    molecule.generate_3d_molecule(opt, plot2D)
    molecule.write_xyz(name + '.xyz')   



# generate input geometries for geometry optimization checking
def generate_xyz_input(smiles_filename, num_smiles, smiles_col_name, id_col_name, cleanup = False, opt = True, plot2D = True):
    """
    generate num_smiles of .xyz files for smiles randomly chosen from csv_file
    :param smiles_filename:
    :param num_smiles:
    :param smiles_col_name: str, the column name of the smiles
    :param id_col_name: str, the column name of the id
    :param cleanup: bool, whether to clean up smiles using chembl_structure_pipeline
    :param opt: bool, whether or not to add H atoms and optimize the MMFF94 force field
    :param plot2D: bool, whether or not to plot 2D structure
    """
    # folder
    save_root_folder = os.path.join(os.getcwd(), 'dataset')
    if not os.path.exists(save_root_folder):
        os.makedirs(save_root_folder)
    
    # generate SMILES_list
    SMILES_list, IDs = rand_smiles_list(smiles_filename, num_smiles, smiles_col_name, id_col_name)
    source = 'smiles_list.txt'
    destination = os.path.join(save_root_folder, source)
    shutil.move(source, destination)
    
    # genertate input .xyz files
    for i, smiles in enumerate(SMILES_list):
        
        if id_col_name:
            name = id_col_name + '_' + str(IDs[i])
        else:
            name = str(IDs[i])
            
        folder = os.path.join(save_root_folder, name)
        if not os.path.exists(folder):
            os.makedirs(folder)
            
        smiles_to_xyz(name, smiles, cleanup, opt, plot2D)
        shutil.move(name + '.xyz', os.path.join(folder, name + '.xyz'))
        if plot2D:
            shutil.move(name + '.png', os.path.join(folder, name + '.png'))
    


if __name__ == '__main__':

    ############################################
    ### Randomly generate SMILES from a file ###
    ############################################
    # smiles_filename = 'TRPML1_ABCFJseries_results.csv'
    # num_smiles = 10
    # smiles_col_name = 'SMILES'
    # id_col_name = 'ID'
    #
    # generate_xyz_input(smiles_filename, num_smiles, smiles_col_name, id_col_name, cleanup = True, opt = True, plot2D = True)
    
    
    ############################################
    ### Generate all SMILES from a file ###
    ############################################
    smiles_filename = 'tests/test_SMILES_to_3D.csv'
    smiles_col_name = 'Cleaned_SMILES'
    id_col_name = 'ID'
    
    process_multiple_smiles(smiles_filename, smiles_col_name, id_col_name, cleanup = True, opt = True, plot2D = True)

