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
from SMILES_to_3D import Molecule


# generate a list of smiles with their IDs
def rand_smiles_list(csv_file, num_smiles):
    """
    randomly choose num_smiles of smiles from input_file_name
    :para csv_file:
    :para num_smiles:
    :return:
    """
    smiles_col_name = 'SMILES'
    id_col_name = 'ID'
    
    data_frame = pd.read_csv(csv_file, header=0)
    sample_num = data_frame.shape[0]
    
    # get smiles and IDs
    if smiles_col_name not in data_frame:
        raise ValueError('SMILES column name "%s" not found in csv file' % smiles_col_name)
    else:
        SMILES_list = data_frame[smiles_col_name].to_numpy()
    
    if id_col_name not in data_frame:
        raise ValueError('ID column name "%s" not found in csv file' % id_col_name)
    else:
        IDs = data_frame[id_col_name].to_numpy()
    
    # shuffle data set
    shuffle_idxs = np.arange(sample_num)
    np.random.shuffle(shuffle_idxs)
    SMILES_list = SMILES_list[shuffle_idxs]
    IDs = IDs[shuffle_idxs]
    
    # get num_smiles of smiles
    SMILES_list = SMILES_list[: num_smiles]
    IDs = IDs[: num_smiles]
    
    with open('smiles_list.txt', 'w') as smiles_list:
        smiles_list.write('ID   original IDs  SMILES \n')
        for i in range(num_smiles):
            outstr = '{}   {}   '.format(i, IDs[i]) + SMILES_list[i] + '\n'
            smiles_list.write(outstr)
    
    return SMILES_list, IDs



# convert a smiles to xyz input
def smiles_to_input(name, smiles):
    """
    convert a smiles to a .xyz file
    :para name: str
    :para smiles: str, smiles string
    :return: None
    """
    assert type(smiles) == str and len(smiles) > 0
    
    molecule = Molecule(name, smiles)
    molecule.smiles_to_3d()
    molecule.output_xyz(name + '.xyz')    



# generate input geometries for geometry optimization checking
def generate_xyz_input(csv_file, num_smiles):
    """
    generate num_smiles of .xyz files for smiles randomly chosen from csv_file
    """
    # folder
    save_root_folder = os.path.join(os.getcwd(), 'dataset')
    if not os.path.exists(save_root_folder):
        os.makedirs(save_root_folder)
    
    # generate SMILES_list
    SMILES_list, IDs = rand_smiles_list(csv_file, num_smiles)
    source = 'smiles_list.txt'
    destination = os.path.join(save_root_folder, source)
    shutil.move(source, destination)
    
    # genertate input .xyz files
    for i in range(num_smiles):
        name = '{}'.format(i)
        folder = os.path.join(save_root_folder, name)
        if not os.path.exists(folder):
            os.makedirs(folder)
            
        smiles_to_input(name, SMILES_list[i])
        shutil.move(name + '.xyz', os.path.join(folder, name + '.xyz'))
    


if __name__ == '__main__':

    csv_file = 'training_phe_antivirus_v3_100nm_7305.csv'
    num_smiles = 10
    generate_xyz_input(csv_file, num_smiles)

    #smiles = 'CC(=O)N[C@@H]1[C@@H](NC(=N)N)C=C(C(=O)O)O[C@H]1[C@H](OC(=O)NCCCCCCCn1cc(CCCCNC(=O)O[C@H]([C@H](O)CO)[C@@H]2OC(C(=O)O)=C[C@H](NC(=N)N)[C@H]2NC(C)=O)nn1)[C@H](O)CO'
    #smiles_to_input('2977', smiles)
