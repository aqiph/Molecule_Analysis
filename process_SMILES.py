#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 10:12:33 2022

@author: guohan

1. process SMILES
2. plot 2D structure
3. clean up SMILES using chembl_structure_pipeline, remove chirality if required
3. clean up disconnected SMILES, i.e., SMILES that contains '.'
4. get molecular features: node features, edge features and topological features

"""

import os, warnings
import numpy as np
import pandas as pd
import pickle
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.Draw import rdMolDraw2D
from chembl_structure_pipeline import *
from chembl_structure_pipeline.checker import *



### process SMILES ###

def process_SMILES(smiles, id, folder, cleanup_SMILES, cleanup_chirality, process_disconnection, process_disconnection_method,
                   addH, plot2D, legend):
    """
    read and process SMILES
    :para smiles: str, SMILES
    :para id: str, the ID of the SMILES
    :para folder: str, the folder path of the output plot
    :para cleanup_SMILES: bool, whether or not to clean up SMILES using chembl_structure_pipeline
    :para cleanup_chirality: bool, whether or not to remove chirality
    :para process_disconnection: bool, whether or not to process disconnected SMILES
    :para process_disconnection_method: str, method for processing other disconnected SMILES,
    if process_disconnection_method == 'keep_longest', keep the longest part in SMILES
    if process_disconnection_method == 'keep_most_atoms', keep the part with the most atoms
    :para plot2D: bool, whether or not plot the SMILES
    :para legend: str or None, legend of plot
    :return: (str, str), SMILES or cleaned SMILES according to 'cleanup' and file name
    """
    # check if SMILES is valid
    if not smiles:
        print(f"Error: Invalid SMILES {smiles}")
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        print(f"Error: Invalid SMILES {smiles}")
        return None
    if not mol:
        print(f"Error: Invalid SMILES {smiles}")
        return None

    # clean up and canonicalize SMILES, remove chirality if required
    if cleanup_SMILES:
        smiles, flag = cleanup_smiles_by_CSP(smiles, cleanup_chirality)
        if not flag:
            smiles = None
    else:
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol)

    # process disconnected SMILES
    if process_disconnection:
        smiles, _ = cleanup_disconnected_smiles(smiles, process_disconnection_method)

    if addH:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        smiles = Chem.MolToSmiles(mol)

    # plot SMILES
    if plot2D:
        if legend == 'ID':
            legend = id
        elif legend == 'SMILES':
            legend = smiles
        else:
            legend = ''

        plot_2d_molecule_from_smiles(smiles, os.path.join(folder, id), legend)

    return smiles, id


def process_single_SMILES(input_file, cleanup_SMILES = True, cleanup_chirality = False,
                          process_disconnection = False, process_disconnection_method = None,
                          addH = False, plot2D = True, legend = None):
    """
    read and process a single SMILES from file
    :para input_file: str, the name of the input file
    :para cleanup_SMILES: bool, whether or not to clean up SMILES using chembl_structure_pipeline
    :para cleanup_chirality: bool, whether or not to remove chirality
    :para process_disconnection: bool, whether or not to process disconnected SMILES
    :para process_disconnection_method: str, method for processing other disconnected SMILES,
    if process_disconnection_method == 'keep_longest', keep the longest part in SMILES
    if process_disconnection_method == 'keep_most_atoms', keep the part with the most atoms
    :para plot2D: bool, whether or not plot the SMILES
    :para legend: str or None, legend of plot
    :return: (str, str), SMILES or cleaned SMILES according to 'cleanup' and file name
    """
    # file name
    folder, basename = os.path.split(os.path.abspath(input_file))
    id = os.path.splitext(basename)[0]
        
    # read SMILES
    with open(input_file, 'r') as f:
        raw = f.readline()
        smiles = raw.strip()

    # process SMILES
    smiles, id = process_SMILES(smiles, id, folder, cleanup_SMILES = cleanup_SMILES, cleanup_chirality = cleanup_chirality,
                                process_disconnection = process_disconnection, process_disconnection_method = process_disconnection_method,
                                addH = addH, plot2D = plot2D, legend = legend)
        
    return smiles, id


def process_multiple_SMILES(input_file, smiles_column_name, id_column_name, cleanup_SMILES = True, cleanup_chirality = False,
                            process_disconnection = False, process_disconnection_method = None,
                            addH = False, plot2D = True, legend = None):
    """
    read and process multiple SMILES from file
    :para input_file: str, the name of the input file
    :para smiles_column_name: str, the column name of the smiles
    :para id_column_name: str, the column name of the id
    :para cleanup_SMILES: bool, whether or not to clean up SMILES using chembl_structure_pipeline
    :para cleanup_chirality: bool, whether or not to remove chirality
    :para process_disconnection: bool, whether or not to process disconnected SMILES
    :para process_disconnection_method: str, method for processing other disconnected SMILES,
    if process_disconnection_method == 'keep_longest', keep the longest part in SMILES
    if process_disconnection_method == 'keep_most_atoms', keep the part with the most atoms
    :para plot2D: bool, whether or not plot the SMILES
    :para legend: str or None, legend of plot
    :return: (str, str), SMILES or cleaned SMILES according to 'cleanup' and file name
    """
    # output name
    folder, basename = os.path.split(os.path.abspath(input_file))
    output_file, fmt = os.path.splitext(basename)

    folder = os.path.join(folder, output_file)
    if not os.path.exists(folder):
        os.makedirs(folder)

    # get DataFrame from input file
    if fmt in {'.xls', '.xlsx'}:
        df = pd.read_excel(input_file)
    elif fmt == '.csv':
        df = pd.read_csv(input_file)
    else:
        print('smiles input file is not readed')
        return

    # get list of id and smiles
    SMILES = df[smiles_column_name].values.tolist()
    num_smiles = len(SMILES)

    if id_column_name is None:
        compound_ID = [i for i in range(num_smiles)]
    else:
        compound_ID = df[id_column_name].values.tolist()

    # generate smiles
    for i, smiles in enumerate(SMILES):
        id = str(compound_ID[i])
        print('Molecule:', id)

        process_SMILES(smiles, id, folder, cleanup_SMILES=cleanup_SMILES, cleanup_chirality=cleanup_chirality,
                       process_disconnection=process_disconnection, process_disconnection_method=process_disconnection_method,
                       addH = addH, plot2D=plot2D, legend=legend)


### plot 2D structure ###

def plot_2d_molecule_from_mol(molecule, output_file_without_ext, legend = '', atomNumber = None):
    """
    plot 2D structure of molecule
    :para molecule: 'mol' object from ridkit
    :para output_file_without_ext: output file name without extension
    :para legend: str, legend of plot
    :para atomNumber: bool or 'mapNumber' or 'atomIndex', add atom map number or add atom indices
    :return: None
    """ 
    # plot molecule       
    d = rdMolDraw2D.MolDraw2DCairo(450, 400)
        
    if atomNumber == 'mapNumber':
        for i, a in enumerate(molecule.GetAtoms()):
            a.SetAtomMapNum(i)
    elif atomNumber == 'atomIndex':
        d.drawOptions().addAtomIndices = True
    
    d.DrawMolecule(molecule, legend = legend)
    
    d.FinishDrawing()
        
    d.WriteDrawingText(output_file_without_ext + '.png')


def plot_2d_molecule_from_smiles(smiles, output_file_without_ext, legend = '', atomNumber = None):
    """
    plot 2D structure of molecule
    :para smiles: str, SMILES of the molecule to plot
    :para output_file_without_ext: output file name without extension
    :para legend: str, legend of plot
    :para atomNumber: bool or 'mapNumber' or 'atomIndex', add atom map number or add atom indices
    :return: None
    """
    # check if SMILES is valid
    if not smiles:
        print(f"Error: Invalid SMILES {smiles}")
        return None
    try:
        molecule = Chem.MolFromSmiles(smiles)
    except Exception:
        print(f"Error: Invalid SMILES {smiles}")
        return None
    if not molecule:
        print(f"Error: Invalid SMILES {smiles}")
        return None
    
    # plot molecule
    d = rdMolDraw2D.MolDraw2DCairo(450, 400)
        
    if atomNumber == 'mapNumber':
        for i, a in enumerate(molecule.GetAtoms()):
            a.SetAtomMapNum(i)
    elif atomNumber == 'atomIndex':
        d.drawOptions().addAtomIndices = True
    
    d.DrawMolecule(molecule, legend = legend)
    
    d.FinishDrawing()
        
    d.WriteDrawingText(output_file_without_ext + '.png')


### clean up SMILES using chembl_structure_pipeline ###

def cleanup_smiles_by_CSP(smiles, cleanup_chirality = False):
    """
    clean up a single SMILES with chembl_structure_pipeline
    :para smiles: str, SMILES
    :para cleanup_chirality: bool, whether or not to remove chirality
    :return: cleaned SMILES by chembl_structure_pipeline, flag to indicate if this SMILES is valid,
    if the SMILES is not valid, return the original SMILES with flag = False
    """
    # True if the input SMILES can be properly converted, False if there is an error
    flag = True

    try:
        # get mol object
        mol = Chem.MolFromSmiles(smiles)
        # standardize mol
        mol_std = standardize_mol(mol)
        # get parent
        mol_par = get_parent_mol(mol_std)[0]
        # canonicalize SMILES, remove chirality if required
        if cleanup_chirality:
            smiles = Chem.MolToSmiles(mol_par, isomericSmiles=False)
            print('Remove chirality')
        else:
            smiles = Chem.MolToSmiles(mol_par)
        mol = Chem.MolFromSmiles(smiles)
        smiles_canonicalized = Chem.MolToSmiles(mol)

    except:
        smiles_canonicalized = smiles
        print(f"Error: Invalid SMILES {smiles}")
        flag = False

    return smiles_canonicalized, flag


### clean up disconnected SMILES, i.e., SMILES that contains '.'

def cleanup_disconnected_smiles(smiles, process_disconnection_method):
    """
    record and process a single disconnected SMILES (containing '.', i.e., polymer, salt, solvent)
    :para smiles: str, SMILES
    :para process_disconnection_method: str, method for processing other disconnected SMILES,
    if process_disconnection_method == 'keep_longest', keep the longest part in SMILES
    if process_disconnection_method == 'keep_most_atoms', keep the part with the most atoms
    """
    if '.' in smiles:
        print('Process disconnected SMILES')

        mol = Chem.MolFromSmiles(smiles)

        # remove simple ions in salts
        ion_list = "[Li,Na,K,Rb,Cs,Mg,Ca,Sr,Ba,Mn,Fe,Co,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,F,Cl,Br,I]"
        saltRemover = SaltRemover(defnData=ion_list)

        mol = saltRemover.StripMol(mol)
        smiles = Chem.MolToSmiles(mol)

        # process other disconnected SMILES
        if '.' in smiles and process_disconnection_method == 'keep_longest':
            smiles = max(smiles.split('.'), key = len)

        elif '.' in smiles and process_disconnection_method == 'keep_most_atoms':
            fragment_list = smiles.split('.')
            fragment_length_list = [len(Chem.MolFromSmiles(fragment).GetAtoms()) for fragment in fragment_list]
            smiles = fragment_list[np.argmax(fragment_length_list)]

        # canonicalize and check SMILES
        mol = Chem.MolFromSmiles(smiles)
        smiles_canonicalized = Chem.MolToSmiles(mol)
        if len(smiles_canonicalized) == 0:
            return None, True
        return smiles_canonicalized, True

    else:
        return smiles, False


### get molecule features ###

def get_feature_from_smiles(smiles, name, addH = False, node_ext_feature=False, plot2D = True, legend = None):
    """
    get molecule features from SMILES
    :para SMILES: str, the SMILES of the input molecule
    :para name: str, the name of the input molecule    
    :para addH: bool, whether or not to add H back to molecule
    :para node_ext_feature: bool, whether or not to add extra features
    :para plot2D: bool, whether or not plot the SMILES
    :para legend: str or None, legend of plot, 
    """
    mol = Chem.MolFromSmiles(smiles)
    
    # get atom_dict
    atom_dict_file = 'utils/atom_dict.pkl'
    atom_dict_file = open(atom_dict_file, "rb")
    atom_dict = pickle.load(atom_dict_file)[0]
    atom_dict_file.close()
    if '<unk>' not in atom_dict:
        atom_dict['<unk>'] = len(atom_dict)
    if False:
        print(atom_dict)
    
    # Canonicalize molecule, add H
    mol = canonicalize_molecule(mol, addH)
    # Get node features
    node_feats = get_node_features(mol, atom_dict, node_ext_feature)
    # Get edge features
    edge_feats = get_edge_features(mol)
    # Get adjacency info
    edge_index = get_adjacency_info(mol)    
    
    # plot smiles
    if plot2D:
        if legend is None:
            legend = name
        elif legend == 'SMILES':
            legend = Chem.MolToSmiles(mol)

        plot_2d_molecule_from_mol(mol, name, legend)
    
    # Get labels info
    return np.transpose(node_feats), np.transpose(edge_feats), edge_index


def canonicalize_molecule(mol, addH = False):
    """
    canonicalize a molecule, add H (optional)
    :param mol: 'mol' object from RDKit
    :param addH: bool, whether or not to add H back to molecule
    :return: molecule canonicalized, 'mol' object from RDKit
    """
    if addH:
        mol = Chem.AddHs(mol)
    order = Chem.CanonicalRankAtoms(mol)
    mol = Chem.RenumberAtoms(mol, order)
    return mol


def get_node_features(mol, atom_dict=None, ext_feature=False):
    all_node_feats = []

    for atom in mol.GetAtoms():
        node_feats = []

        # Feature 1: Atomic number from dict
        if atom_dict:
            if atom.GetSymbol() not in atom_dict:
                node_feats.append(atom_dict['<unk>'])
                warnings.warn(f'{atom.GetSymbol()} is not in atom_dict')
            else:
                node_feats.append(atom_dict[atom.GetSymbol()])
        else:
            node_feats.append(atom.GetAtomicNum())

        if ext_feature:
            # Feature 2: Atom degree
            node_feats.append(atom.GetDegree())
            # Feature 3: Formal charge
            node_feats.append(atom.GetFormalCharge())
            # Feature 4: Hybridization
            node_feats.append(atom.GetHybridization())
            # Feature 5: Aromaticity
            node_feats.append(atom.GetIsAromatic())
            # Feature 6: Total Num Hs
            node_feats.append(atom.GetTotalNumHs())
            # Feature 7: Radical Electrons
            node_feats.append(atom.GetNumRadicalElectrons())
            # Feature 8: In Ring
            node_feats.append(atom.IsInRing())
            # Feature 9: Chirality
            node_feats.append(atom.GetChiralTag())

        # Append node features to matrix
        all_node_feats.append(node_feats)

    all_node_feats = np.asarray(all_node_feats)
    return all_node_feats


def get_edge_features(mol):
    """
    This will return a matrix / 2d array of the shape
    [Number of edges, Edge Feature size]
    :param mol: 'mol' object from RDKit
    """
    all_edge_feats = []

    for bond in mol.GetBonds():
        edge_feats = []
        # Feature 1: Bond type (as double)
        edge_feats.append(bond.GetBondTypeAsDouble())
        # Feature 2: Rings
        edge_feats.append(bond.IsInRing())
        # Append node features to matrix (twice, per direction)
        all_edge_feats += [edge_feats, edge_feats]

    all_edge_feats = np.asarray(all_edge_feats)
    return all_edge_feats


def get_adjacency_info(mol):
    """
    We could also use rdmolops.GetAdjacencyMatrix(mol)
    but we want to be sure that the order of the indices
    matches the order of the edge features
    :param mol: 'mol' object from RDKit
    """
    edge_indices = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edge_indices += [[i, j], [j, i]]

    edge_indices = np.asarray(edge_indices)
    edge_indices = np.transpose(edge_indices)
    return edge_indices



####

if __name__ == '__main__':

    multiple = False

    if not multiple:
        ############################################
        ### process one smile from an input file ###
        ############################################
        input_file = 'process_SMILES/tests/test.txt'
        smiles, id = process_single_SMILES(input_file, cleanup_SMILES = True, cleanup_chirality = False,
                                           process_disconnection = False, process_disconnection_method = 'keep_most_atoms',
                                           addH = False, plot2D = True, legend = 'SMILES')

        features = False
        if features:
            print('#########')
            id += '_get_feature'
            node_feats, edge_feats, edge_index = get_feature_from_smiles(smiles, id, addH = False, node_ext_feature = True, plot2D = True, legend = 'SMILES')
        
            print('#########')
            print('Node features')
            print(node_feats)
            print('#########')
            print('Edge features')
            print(edge_feats)
            print('#########')
            print('Edge index')
            print(edge_index)
    
    else:
        ##########################################
        ### process multiple smiles from file ####
        ##########################################
        input_file = 'process_SMILES/tests/test.csv'
        smiles_column_name = 'SMILES'
        id_column_name = 'ID'

        process_multiple_SMILES(input_file, smiles_column_name, id_column_name, cleanup_SMILES=False, cleanup_chirality=False,
                                process_disconnection=False, process_disconnection_method='keep_most_atoms',
                                addH = False, plot2D=True, legend='SMILES')


