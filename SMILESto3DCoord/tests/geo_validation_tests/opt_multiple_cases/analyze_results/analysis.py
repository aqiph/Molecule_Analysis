#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 16:34:17 2021

@author: guohan
"""
import os, time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_family('serif')
font.set_name('Calibri')
font.set_size(13)

import ase
from ase.io import read, write
from ase import Atoms



class Molecule(object):
    
    
    def __init__(self, in_filename = None):
        
        assert isinstance(in_filename, str)           
        atoms = read(in_filename, format = 'xyz')
        assert isinstance (atoms, Atoms)
        
        self.atoms = atoms
        self.positions = atoms.get_positions()
        self.atom_list = atoms.get_chemical_symbols()
        self.num_atoms = len(self.atom_list)
        
        self.distance_matrix = []
        
        for i in range(self.num_atoms):
            distance_list = []
            atom1 = self.positions[i]
            for j in range(self.num_atoms):
                if i == j:
                    distance_list.append(1000.0)
                else:
                    atom2 = self.positions[j]
                    distance_list.append(np.linalg.norm(atom2 - atom1))
            self.distance_matrix.append(distance_list)
        self.distance_matrix = np.array(self.distance_matrix)
        
        assert len(self.distance_matrix) == self.num_atoms
        assert len(self.distance_matrix[0]) == self.num_atoms
    
    
    def connectivity_mask(self, cutoff_radius):
        """
        compute the connectivity mask accordinging to 'cutoff_radius'
        :para cutoff_radius: float, the cutt radius 
        :return mask: matrix of bool, True if the bondlength <= cutoff_radius, False otherwise
        """
        mask = [self.distance_matrix <= cutoff_radius]
        return mask
    
    
    def bondlength_change_rmsd(self, other, cutoff_radius):
        """
        compute the rmsd of the bondlength changes accordinging to 'cutoff_radius'
        :para other: a Molecule object
        :para cutoff_radius: float, the cutt radius 
        :return:
        """
        if not isinstance(other, Molecule):
            print('Error: input object is not an instance of Molecule')
        
        mask = self.connectivity_mask(cutoff_radius)
        num_bonds = np.sum(mask)
        rmsd = np.sum(((self.distance_matrix - other.distance_matrix) * mask) ** 2.0) / (2 * num_bonds)
        
        return rmsd
        
    
    def bondlength_change(self, other, cutoff_radius, large_discrepancy = 0.5):
        """
        compute the bondlength change accordinging to 'cutoff_radius'
        :para other: a Molecule object
        :para cutoff_radius: float, the cutt radius 
        :return:
        """
        if not isinstance(other, Molecule):
            print('Error: input object is not an instance of Molecule')
        
        num_bonds = 0
        bondlength_changes = []
        bond_types = []
        bondlength_change_dict = {}
        
        bondlength_change_matrix = self.distance_matrix - other.distance_matrix
        
        for i in range(self.num_atoms):
            atom1_type = self.atom_list[i]
            
            for j in range(i + 1, self.num_atoms):
                if self.distance_matrix[i][j] <= cutoff_radius:
                    num_bonds += 1
                    bondlength_change = bondlength_change_matrix[i][j]
                    bondlength_changes.append(bondlength_change)  
                    if bondlength_change > large_discrepancy or bondlength_change < -large_discrepancy:
                        print('The change in Atom{}-Atom{} is large {}'.format(i, j, bondlength_change))
                    
                    atom2_type = self.atom_list[j]
                    bond_type = sorted([atom1_type, atom2_type])
                    bond_type = bond_type[0] + '-' + bond_type[1]
                    bond_types.append(bond_type)
                    
                    if bond_type not in bondlength_change_dict:
                        bondlength_change_dict[bond_type] = []
                    bondlength_change_dict[bond_type].append(bondlength_change)
        
#        print(bondlength_changes)
#        print(bond_types)
#        print(bondlength_change_dict)
        
        return bondlength_changes, bond_types, bondlength_change_dict
    
    
def plot_bondlength_change_dist(cutoff_radius, bond_type = None, large_discrepancy = 0.5):
    """
    plot distribution of bondlength change
    :para cutoff_radius:
    :para bodn_type: str, 
    """
    bondlength_changes = []
    
    cwd = os.getcwd()
    
    # compute bond length change
    for i in range(10):
        print('Molecule {}:'.format(i))
        in_filename = os.path.join(cwd, 'results', str(i), str(i) + '.xyz')
        molecule = Molecule(in_filename = in_filename)
        in_filename_opt = os.path.join(cwd, 'results', str(i), str(i) + '_opt.xyz')
        molecule_opt = Molecule(in_filename = in_filename_opt)
        
        bondlength_changes_all, bond_types, bondlength_change_dict = molecule.bondlength_change(molecule_opt, cutoff_radius)
    
        if bond_type == 'all' or bond_type is None:
            bondlength_changes.extend(bondlength_changes_all)
            title = 'all bonds'
        else:
            if bond_type not in bondlength_change_dict:
                print('Bond type not found!')
                return
            bondlength_changes.extend(bondlength_change_dict[bond_type])
            title = bond_type
    
    # number of bonds with large discrepancy
    num_bonds_all = len(bondlength_changes)
    num_extreme_bonds_all = np.sum(np.array(bondlength_changes) > large_discrepancy) + np.sum(np.array(bondlength_changes) < -large_discrepancy)
    
    # statistics
    ave = np.average(np.array(bondlength_changes))
    std = np.std(np.array(bondlength_changes))
    
    # plot
    plt.figure(1)
    low_range = -large_discrepancy
    high_range = large_discrepancy
    stepsize = (high_range - low_range) / 20.0
    bins = np.array([low_range + i * stepsize for i in range(20)])
    plt.hist(bondlength_changes, bins = bins)
    
    plt.xlabel('Bond length change (' + r'$\mathrm{\AA}$' + ')', fontproperties = font)
    plt.ylabel('Number of bonds', fontproperties = font)
    plt.xticks(fontproperties = font)
    plt.yticks(fontproperties = font)
    plt.xlim(low_range, high_range)
#    plt.legend(frameon = False, fontsize = 15)
    plt.title('Bond type: ' + title + ' ({:3.2f} %)'.format(100.0 - 100.0*num_extreme_bonds_all/num_bonds_all), fontproperties = font)
    
    plt.text(-0.4, 160, 'ave: {:2.6f}\nstd: {:2.6f}'.format(ave, std))
    
    plt.savefig('bondlength_change_dist.png', dpi = 300, bbox_inches = 'tight')
    plt.show()
    plt.close()
        
                    
        
        


if __name__ == '__main__':
    
    cutoff_radius = 7.0
    
    plot_bondlength_change_dist(cutoff_radius, bond_type = 'H-O')
    
    
    
    
    