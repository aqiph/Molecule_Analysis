#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 10:47:06 2021

@author: guohan



"""

import os
import shutil
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from process_SMILES import plot_2d_molecule_from_mol, cleanup_smiles_by_CSP



def process_multiple_smiles(input_file, smiles_col_name, id_col_name = None, cleanup_SMILES = False, opt = True, plot2D = True):
    """
    read smiles from excel, generate .xyz file of each smiles in the input_file
    :para input_file: str, name of input file
    :para smiles_col_name: str, the column name of the smiles
    :para id_col_name: str, the column name of the id
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
    
    # get list of id and smiles
    SMILES = df[smiles_col_name].values.tolist()
    num_smiles = len(SMILES)
    
    if id_col_name is None:
        compound_ID = [i for i in range(num_smiles)]
    else:
        compound_ID = df[id_col_name].values.tolist()
    
    # generate smiles
    for i, smiles in enumerate(SMILES):
        print('Molecule:', i)
        
        if id_col_name:
            name = id_col_name + '_' + str(compound_ID[i])
        else:
            name = str(compound_ID[i])
        
        # cleanup SMILES
        if cleanup_SMILES:
            smiles, flag = cleanup_smiles_by_CSP(smiles, cleanup_chirality = False)
            if not flag:
                smiles = None
        
        # generate Molecule object and 3D structure
        molecule = Molecule(name, smiles)    
        molecule.generate_3d_molecule(opt, plot2D)
        if plot2D:
            shutil.move(name + '.png', os.path.join(folder, name + '.png'))
        
        molecule.write_xyz(name + '.xyz')        
        shutil.move(name + '.xyz', os.path.join(folder, name + '.xyz'))
        
        with open(name + '.txt', 'w') as f:
            f.write(smiles)        
        shutil.move(name + '.txt', os.path.join(folder, name + '.txt'))


class Molecule(object):
    
    def __init__(self, name = None, smiles = None):  
        self.folder = os.getcwd()
        self.name = name
        self.smiles = smiles
        self.xyz = []
        self.atoms = []
        self.num_atoms = 0
        
        assert type(smiles) == str or smiles is None, 'SMILES is not valid'    
        
    
    def set_smiles(self, smiles):        
        assert type(smiles) == str, 'SMILES is not valid'        
        self.smiles = smiles
        
    
    def get_smiles(self):       
        return self.smiles
        
    
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
    
    
    def generate_3d_molecule(self, opt = True, plot2D = True):
        """
        Generate .xyz file of the molecule. If opt == True, optimize molecule with the MMFF94 force field
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
            plot_2d_molecule_from_mol(molecule, output_file_without_ext = os.path.join(self.folder, self.name), legend = self.name, atomNumber = None)
        
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
    
    
    def read_smiles(self, input_file, cleanup_SMILES = False):
        """
        read SMILES from file
        :para input_file: str, file name
        :para cleanup_SMILES: bool, whether or not to clean up SMILES using chembl_structure_pipeline
        :return: str, SMILES or cleaned SMILES according to 'cleanup'
        """  
        # file name        
        folder, basename = os.path.split(os.path.abspath(input_file))
        output_file, fmt = os.path.splitext(basename)
        self.folder = folder
        self.name = output_file
        
        # read smiles 
        with open(input_file, 'r') as f:
            raw = f.readline()
            smiles = raw.strip()
        
        # cleanup smiles
        if cleanup_SMILES:
            cleaned_SMILES, flag = cleanup_smiles_by_CSP(smiles, cleanup_chirality = False)
            if not flag:
                cleaned_SMILES = None
                
            self.smiles = cleaned_SMILES
            return cleaned_SMILES
        
        self.smiles = smiles    
        return smiles
        
    
    def write_xyz(self, output_file = None):
        """
        Output .xyz file
        :para output_file: str, output file name
        """        
        if len(self.xyz) <= 0:
            print('xyz is empty')
            return
        
        if output_file is None:
            output_file = os.path.join(self.folder, self.name + '.xyz')
        
        f = open(output_file, 'w')
        f.write(str(self.num_atoms) + '\n')
        f.write(self.smiles + '\n')
        
        for n in range(self.num_atoms):
            atom_xyz = self.xyz[n]
            out_str = self.atoms[n] + '   {:4.4f}   {:4.4f}   {:4.4f} \n'.format(atom_xyz[0], atom_xyz[1], atom_xyz[2])
            f.write(out_str)
            


if __name__ == '__main__':
    
    multiple = True
    
    if not multiple:
        ############################################
        ### process one smile from an input file ###
        ############################################
        input_file = 'SMILES_to_3D/tests/test.txt'

        molecule = Molecule()   
        molecule.read_smiles(input_file, cleanup_SMILES = True)
        molecule.generate_3d_molecule(opt = True, plot2D = True)
        molecule.write_xyz()
    
    else:    
        ##########################################
        ### process multiple smiles from file ####
        ##########################################
        input_file = 'SMILES_to_3D/tests/test_SMILES_to_3D.csv'
        smiles_col_name = 'SMILES'
        id_col_name = 'ID'

        process_multiple_smiles(input_file, smiles_col_name, id_col_name, cleanup_SMILES = True, opt = True, plot2D = True)
    
    
