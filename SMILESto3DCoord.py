#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 10:47:06 2021

@author: guohan

1. convert single SMILES into 3D structure
2. convert multiple SMILES into 3D structure

"""

import os
import shutil
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from process_SMILES import plot_2d_molecule_from_mol, cleanup_smiles_by_CSP



### convert SMILES to 3D coordinates ###

def convert_single_SMILES_to_Coord(input_file, cleanup_SMILES = False, opt = True, plot2D = True):
    """
    read a single SMILES from excel, generate .xyz file of each smiles in the input_file
    :para input_file: str, the name of the input file
    :para cleanup_SMILES: bool, whether or not to clean up SMILES using chembl_structure_pipeline
    :para opt: bool, whether or not to add H atoms and optimize the MMFF94 force field
    :para plot2D: bool, whether or not to plot 2D structure
    """
    # output name
    folder, basename = os.path.split(os.path.abspath(input_file))
    output_file, fmt = os.path.splitext(basename)

    # get SMILES and ID
    id = output_file

    with open(input_file, 'r') as f:
        raw = f.readline()
        smiles = raw.strip()

    # cleanup SMILES
    if cleanup_SMILES:
        smiles, flag = cleanup_smiles_by_CSP(smiles, cleanup_chirality=False)
        if not flag:
            smiles = None

    # generate 3D coordinates
    molecule = Molecule(id, smiles)
    molecule.generate_3d_molecule(folder, opt, plot2D)
    molecule.write_xyz(folder)


def convert_multiple_SMILES_to_Coord(input_file, smiles_column_name, id_column_name = None, cleanup_SMILES = False, opt = True, plot2D = True):
    """
    read multiple SMILES from excel, generate .xyz file of each smiles in the input_file
    :para input_file: str, the name of the input file
    :para smiles_column_name: str, the column name of the smiles
    :para id_column_name: str, the column name of the id
    :para cleanup_SMILES: bool, whether or not to clean up SMILES using chembl_structure_pipeline
    :para opt: bool, whether or not to add H atoms and optimize the MMFF94 force field
    :para plot2D: bool, whether or not to plot 2D structure
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
    
    # get SMILES and IDs
    SMILES = df[smiles_column_name].values.tolist()
    num_smiles = len(SMILES)
    
    if id_column_name is None:
        compound_ID = [i for i in range(num_smiles)]
    else:
        compound_ID = df[id_column_name].values.tolist()

    # loop
    for i, smiles in enumerate(SMILES):
        id = str(compound_ID[i])
        print('Molecule:', id)
        
        # cleanup SMILES
        if cleanup_SMILES:
            smiles, flag = cleanup_smiles_by_CSP(smiles, cleanup_chirality = False)
            if not flag:
                smiles = None
        
        # generate 3D coordinates
        molecule = Molecule(id, smiles)
        molecule.generate_3d_molecule(folder, opt, plot2D)
        molecule.write_xyz(folder)


class Molecule(object):
    
    def __init__(self, id = None, smiles = None):
        self.id = id
        self.smiles = smiles
        self.xyz = []
        self.atoms = []
        self.num_atoms = 0
        
        assert type(smiles) == str or smiles is None, 'SMILES is not valid'

    
    def add_H(self, smiles):
        """
        canonicalize smiles, convert to mol object
        add H -> canonicalize molecule with modified order
        :param smiles: smiles to be converted
        :return: molecule canonicalized, 'mol' object from rdkit
        """        
        # convert SMILES to mol object in rdkit
        mol = Chem.MolFromSmiles(smiles)
        
        # add H
        mol_with_H = Chem.AddHs(mol)
        
        # canonicalize molecule with modified order
        neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol_with_H))])))[1]        
        mol_with_H = Chem.RenumberAtoms(mol_with_H, neworder)
        
        return mol_with_H
    
    
    def generate_3d_molecule(self, folder, opt = True, plot2D = True):
        """
        Generate .xyz file of the molecule. If opt == True, optimize molecule with the MMFF94 force field
        :para folder: str, the name of the folder
        :para opt: bool, whether or not to add H atoms and optimize the MMFF94 force field
        :para plot2D: bool, whether or not to plot 2D structure
        :return: None
        """        
        if self.smiles is None:
            print('Need to input smiles')
            return
        
        if opt:        
            # add H, canonicalize molecule with modified order, update self.smiles
            molecule = self.add_H(self.smiles)                
            self.smiles = Chem.MolToSmiles(molecule, isomericSmiles=True, kekuleSmiles=True,
                                           canonical=True, allBondsExplicit=True, allHsExplicit=True)
        
            # embed molecule using ETKDG method and optimize molecule with the MMFF94 force field
            AllChem.EmbedMolecule(molecule)
            AllChem.MMFFOptimizeMolecule(molecule)
            
            self.smiles = Chem.MolToSmiles(molecule, isomericSmiles=True, kekuleSmiles=True,
                                           canonical=True, allBondsExplicit=True, allHsExplicit=True)
        
        else:
            # get mol object, update self.smiles
            molecule = Chem.MolFromSmiles(self.smiles)
            self.smiles = Chem.MolToSmiles(molecule, isomericSmiles=True, kekuleSmiles=True,
                                           canonical=True, allBondsExplicit=True)
            molecule = Chem.MolFromSmiles(self.smiles)
        
        # plot 2d structure if needed
        if plot2D:
            plot_2d_molecule_from_mol(molecule, output_file_without_ext = os.path.join(folder, self.id), legend = self.id, atomNumber = None)
        
        # get xyz of molecule
        mol_out = Chem.MolToMolBlock(molecule)
        
        # generate to xyz file
        mol_out = mol_out.split('\n')        
        self.num_atoms = molecule.GetNumAtoms()
        
        if False:
            print('number of atoms', self.num_atoms)
        
        for n in range(self.num_atoms):
            line = mol_out[n + 4].split()
            assert len(line) >= 4
            
            atom_xyz = [float(i) for i in line[:3]]
            atom_type = line[3]
            
            self.xyz.append(atom_xyz)
            self.atoms.append(atom_type)

    
    def write_xyz(self, folder):
        """
        output .xyz file
        :para folder: str, the name of the folder
        """        
        if len(self.xyz) <= 0:
            print('xyz is empty')
            return
        
        output_file = os.path.join(folder, self.id + '.xyz')

        f = open(output_file, 'w')
        f.write(str(self.num_atoms) + '\n')
        f.write(self.smiles + '\n')
        
        for n in range(self.num_atoms):
            atom_xyz = self.xyz[n]
            out_str = self.atoms[n] + '   {:4.4f}   {:4.4f}   {:4.4f} \n'.format(atom_xyz[0], atom_xyz[1], atom_xyz[2])
            f.write(out_str)
            


if __name__ == '__main__':
    
    multiple = False
    
    if not multiple:
        ############################################
        ### process one smile from an input file ###
        ############################################
        input_file = 'SMILESto3DCoord/tests/test_convert_single_SMILES_to_Coord.txt'

        convert_single_SMILES_to_Coord(input_file, cleanup_SMILES = False, opt = True, plot2D = True)
    
    else:    
        ##########################################
        ### process multiple smiles from file ####
        ##########################################
        input_file = 'SMILESto3DCoord/tests/test_convert_multiple_SMILES_to_Coord.csv'
        smiles_column_name = 'SMILES'
        id_column_name = 'ID'

        convert_multiple_SMILES_to_Coord(input_file, smiles_column_name, id_column_name, cleanup_SMILES = False, opt = True, plot2D = True)
    
    
