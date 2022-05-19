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
from rdkit.Chem.Draw import rdMolDraw2D

from chembl_structure_pipeline import *
from chembl_structure_pipeline.checker import *


def preprocess_smiles(smiles):
    """
    helper function for cleanning up smiles using cleanup_smiles_GChem
    :para smiles: str smiles
    :return: cleaned smiles by chembl_structure_pipeline
    """
    flag = 1
    try:
        mol = Chem.MolFromSmiles(smiles)
        # Standardize mol
        mol_std = standardize_mol(mol)
        standardized_smiles = Chem.MolToSmiles(mol_std)
        # get parent
        mol_par = get_parent_mol(mol_std)[0]
        # mol to SMILES
        canonical_smiles = Chem.MolToSmiles(mol_par)
        mol_canonical = Chem.MolFromSmiles(canonical_smiles)
        canonical_smiles = Chem.MolToSmiles(mol_canonical)
        
    except:
        standardized_smiles = smiles
        canonical_smiles = smiles
        print(f"Error for SMILES: {smiles}")
        flag = 0
        
    return canonical_smiles, flag


def plot_2d_molecule(molecule, name, legend = '', atomNumber = None):
    """
    plot 2D structure of molecule
    :para molecule: 'mol' object from ridkit
    :para name: output file name without extension
    :para legend: str, legend of plot
    :para atomNumber: bool or 'mapNumber' or 'atomIndex', add atom map number or add atom indices
    :return: None
    """        
    d = rdMolDraw2D.MolDraw2DCairo(450, 400)
        
    if atomNumber == 'mapNumber':
        for i, a in enumerate(molecule.GetAtoms()):
            a.SetAtomMapNum(i)
    elif atomNumber == 'atomIndex':
        d.drawOptions().addAtomIndices = True
    
    d.DrawMolecule(molecule, legend = legend)
    
    d.FinishDrawing()
        
    d.WriteDrawingText(name + '.png')


def process_multiple_smiles(smiles_filename, smiles_name, id_name = None):
    """
    read smiles from excel, generate .xyz file of each smiles in the smiles_filename
    :para smiles_filename: str, name of input file
    :para smiles_name: str, the column name of the smiles
    :para id_name: str, the column name of the id
    """
    # file name
    _, basename = os.path.split(os.path.abspath(smiles_filename))
    basename_without_ext  = os.path.splitext(basename)[0]
    
    output_folder = os.path.join(os.getcwd(), basename_without_ext)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # get DataFrame from excel
    df = pd.read_excel(smiles_filename)
    
    # get list of id and smiles
    SMILES = df[smiles_name].values.tolist()
    num_smiles = len(SMILES)
    
    if id_name is None:
        compound_ID = [i for i in range(num_smiles)]
    else:
        compound_ID = df[id_name].values.tolist()
    
    # generate smiles
    for i, smiles in enumerate(SMILES):
        print('Molecule:', i)
        
        name = str(compound_ID[i])
        molecule = Molecule(name, smiles)    
        molecule.generate_3d_molecule(addH = True, plot2D = False)
        molecule.write_xyz(name + '.xyz')
        
        shutil.move(name + '.xyz', os.path.join(output_folder, name + '.xyz'))
        
        with open(name + '.txt', 'w') as f:
            f.write(smiles)
        
        shutil.move(name + '.txt', os.path.join(output_folder, name + '.txt'))


class Molecule(object):
    
    def __init__(self, name = None, smiles = None):        
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
        
    
    def add_H_v0(self, smiles):
        """
        canonicalize smiles, convert to mol object
        add H -> canonicalize molecule
        :param smiles: smiles to be converted
        :return: molecule canonicalized, 'mol' object from rdkit
        """       
        # convert SMILES to mol object in rdkit
        mol = Chem.MolFromSmiles(smiles)
        
        # add H
        mol_with_H = Chem.AddHs(mol)        
        self.plot_2d_molecule(mol_with_H, name = 'fig_before_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        # canonicalize molecule with order
        order = Chem.CanonicalRankAtoms(mol_with_H)
        mol_with_H = Chem.RenumberAtoms(mol_with_H, order)        
        self.plot_2d_molecule(mol_with_H, name = 'fig_after_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        return mol_with_H
    
    
    def add_H_v1(self, smiles):
        """
        canonicalize smiles, convert to mol object
        canonicalize molecule -> add H
        :param smiles: smiles to be converted
        :return: molecule canonicalized, 'mol' object from rdkit
        """        
        # convert SMILES to mol object in rdkit
        mol = Chem.MolFromSmiles(smiles)        
        self.plot_2d_molecule(mol, name = 'fig_before_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        # canonicalize molecule with order
        order = Chem.CanonicalRankAtoms(mol)
        mol = Chem.RenumberAtoms(mol, order)
        
        # add H
        mol_with_H = Chem.AddHs(mol)        
        self.plot_2d_molecule(mol_with_H, name = 'fig_after_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        return mol_with_H
    
    
    def add_H_v2(self, smiles):
        """
        canonicalize smiles, convert to mol object
        add H -> canonicalize molecule -> convert to smiles, convert to 'mol' object
        :param smiles: smiles to be converted
        :return: molecule canonicalized, 'mol' object from rdkit
        """        
        # convert SMILES to mol object in rdkit
        mol = Chem.MolFromSmiles(smiles)
        
        # add H
        mol_with_H = Chem.AddHs(mol)        
        self.plot_2d_molecule(mol_with_H, name = 'fig_before_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        # canonicalize molecule with order
        order = Chem.CanonicalRankAtoms(mol_with_H)
        mol_with_H = Chem.RenumberAtoms(mol_with_H, order)        
        self.plot_2d_molecule(mol_with_H, name = 'fig_afterReorder_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        # convert to smiles, then convert back to 'mol' object
        smiles = Chem.MolToSmiles(mol_with_H, isomericSmiles=True, kekuleSmiles=True,
                                  canonical=True, allBondsExplicit=True, allHsExplicit=True)
        molecule = Chem.MolFromSmiles(smiles)
        molecule = Chem.AddHs(molecule)        
        self.plot_2d_molecule(molecule, name = 'fig_fin_' + self.name, legend = self.name, atomNumber = 'atomIndex')        
        
        return molecule
    
    
    def add_H_v3(self, smiles):
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
        self.plot_2d_molecule(mol_with_H, name = 'fig_before_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        # canonicalize molecule with modified order
        neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol_with_H))])))[1]     
        mol_with_H = Chem.RenumberAtoms(mol_with_H, neworder)        
        self.plot_2d_molecule(mol_with_H, name = 'fig_after_' + self.name, legend = self.name, atomNumber = 'atomIndex')
        
        return mol_with_H
    
    
    def generate_3d_molecule(self, addH = True, plot2D = True):
        """
        Generate .xyz file of the molecule. If addH == True, optimize molecule with the MMFF94 force field
        :para addH: bool, whether or not to add H atoms and optimize the MMFF94 force field
        :para plot2D: bool, whether or not to plot 2D structure
        :return: None
        """
        if self.smiles is None:
            print('Need to input smiles')
            return
        
        if addH:        
            # add H, canonicalize molecule with order, update self.smiles
            molecule = self.add_H_v3(self.smiles)                
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
            self.plot_2d_molecule(molecule, name = self.name, legend = self.name, atomNumber = 'atomIndex')
            
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
    
    
    def read_smiles(self, input_filename, cleanup = False):
        """
        read smiles from file
        :para input_filename: str, file name
        :para cleanup: bool, whether to clean up smiles using chembl_structure_pipeline
        :return: str, smiles or cleaned smiles according to 'cleanup'
        """  
        # file name
        _, basename = os.path.split(os.path.abspath(input_filename))
        basename_without_ext  = os.path.splitext(basename)[0]
        self.name = basename_without_ext
        
        # read smiles 
        with open(input_filename, 'r') as f:
            raw = f.readline()
            smiles = raw.strip()
        
            # cleanup smiles
            if cleanup:
                smiles_cleanup, flag = preprocess_smiles(smiles)
                
                if False:
                    print(smiles_cleanup)
                
                if flag == 0:
                    print('SMILES is not correct')
                    return smiles_cleanup
                
                self.smiles = smiles_cleanup
                return smiles_cleanup
        
        self.smiles = smiles    
        return smiles
        
    
    def write_xyz(self, output_filename = None):
        """
        Output .xyz file
        :para output_filename: str, output file name
        """        
        if len(self.xyz) <= 0:
            print('xyz is empty')
            return
        
        if output_filename is None:
            output_filename = os.path.join(os.getcwd(), self.name + '.xyz')
        
        f = open(output_filename, 'w')
        f.write(str(self.num_atoms) + '\n')
        f.write(self.smiles + '\n')
        
        for n in range(self.num_atoms):
            atom_xyz = self.xyz[n]
            out_str = self.atoms[n] + '   {:4.4f}   {:4.4f}   {:4.4f} \n'.format(atom_xyz[0], atom_xyz[1], atom_xyz[2])
            f.write(out_str)
            


if __name__ == '__main__':
    
    # process one smile from an input file
    input_filename = 'test.txt'

    molecule = Molecule()   
    molecule.read_smiles(input_filename, cleanup = True)
    molecule.generate_3d_molecule(addH = False, plot2D = True)
    molecule.write_xyz()
    
    # process multiple smiles from excel
    # smiles_filename = 'pi-pi stacking.xlsx'
    
    # process_multiple_smiles(smiles_filename, smiles_name = 'SMILES', id_name = 'Compound_ID')
