#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 15:25:47 2021

@author: guohan

1. fingerprint
2. scaffold

"""

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem.Scaffolds import MurckoScaffold
# from mhfp.encoder import MHFPEncoder



def get_fingerprint(smiles, fp_method = 'ecfp4'):
    """
    compute the fingerprint of the input smiles
    :para smiles: str, input SMILES
    :para fp_method: str, method to compute fingerprint (topology: topological fingerprint; 
                                                         maccs: MACCS fingerprint;
                                                         atomPairs: atom pairs;
                                                         torsions: topological torsions;
                                                         ecfp4: Morgan ECFP4 (connectivity) fingerprint;
                                                         fcfp4: Morgan FCFP4 (feature-based) fingerprint)
                                                         mhfp6: MinHash fingerprint
    :return: fingerprint object, list, or None
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
    
    # calculate fingerprint
    fp = None
    fp_method = fp_method.lower()
    
    if fp_method == 'topology':
        fp = Chem.RDKFingerprint(mol)
    elif fp_method == 'maccs':
        fp = MACCSkeys.GenMACCSKeys(mol)
    elif fp_method == 'atompairs':
        fp = Pairs.GetAtomPairFingerprint(mol)
    elif fp_method == 'torsions':
        fp = Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)
    elif fp_method == 'ecfp4':
        fp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, )
    elif fp_method == 'fcfp4':
        fp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True)
    # elif fp_method == 'mhfp6':
    #     mhfp_encoder = MHFPEncoder(n_permutations =  2048, seed = 42)
    #     fp = mhfp_encoder.encode(smiles)
    else:
        raise NotImplementedError(f"Fingerprint method {fp_method} not implemented")
    
    return fp


def cal_fingerprint_distance(fp_1, fp_2, fp_method = 'ecfp4'):
    """
    use Tanimoto/Jaccard algorithm to compute distance between fp_1 and fp2: 
    similarity = #(common feature)/len(fp), distance = 1 - similarity
    :para fp_1: fingerprint for molecule_1
    :para fp_2: fingerprint for molecule_2
    :para fp_method: str, method to compute fingerprint (topology: topological fingerprint;
                                                         maccs: MACCS fingerprint;
                                                         atomPairs: atom pairs;
                                                         torsions: topological torsions;
                                                         ecfp4: Morgan ECFP4 (connectivity) fingerprint;
                                                         fcfp4: Morgan FCFP4 (feature-based) fingerprint)
                                                         mhfp6: MinHash fingerprint
    :return: int or None
    """
    if fp_1 is None or fp_2 is None:
        return 1.0
    
    fp_method = fp_method.lower()
    
    if fp_method in {'topology', 'maccs', 'atompairs', 'torsions', 'ecfp4', 'fcfp4'}:
        sim = DataStructs.TanimotoSimilarity(fp_1, fp_2)
        dist = 1.0 - sim
    
    # elif fp_method == 'mhfp6':
    #     dist = MHFPEncoder.distance(fp_1, fp_2)
    
    else:
        raise NotImplementedError(f"Fingerprint method {fp_method} not implemented")
    
    return dist    


def get_scaffold(smiles, include_chirality = False, generic = False):
    """
    compute the Bemis-Murcko scaffold for a SMILES string
    :para smiles: str, input SMILES
    :para include_chirality: bool, whether or not using chirality
    :para generic: bool, whether or not make the scaffold generic
    (i.e., change all heavy atoms to the same one, ignore bond order)
    :return: str or None, SMILES of the BM scaffold
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
    
    # calculate scaffold
    if generic == True:
        scaffold_mol = MurckoScaffold.MakeScaffoldGeneric(mol=mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold_mol)
    else:
        scaffold_smiles = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=include_chirality)
    
    return scaffold_smiles



if __name__ == '__main__':
    
    ### test get_fingerprint
    # smiles = 'Cc1cc(Oc2nccc(CCC)c2)ccc1'
    # smiles = ''
    # plot_2d_molecule_from_smiles(smiles, output_file_without_ext = 'original', legend = smiles, atomNumber = None)
    # print(get_fingerprint(smiles, fp_method = 'MHFP6'))
    
    ### test cal_fingerprint_distance
    smiles_1 = 'CCOC1=C(C=C(C=C1)S(=O)(=O)N(C)C)C2=NC(=O)C3=C(N2)C(=NN3C)C(C)(C)C'
    smiles_2 = 'CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C'
    # smiles_1 = 'CCO'
    # smiles_2 = 'h'

    for fp_method in ['topology', 'maccs', 'atompairs', 'torsions', 'ecfp4', 'fcfp4']:
        fp_1 = get_fingerprint(smiles_1, fp_method = fp_method)
        fp_2 = get_fingerprint(smiles_2, fp_method = fp_method)
        print(fp_method, cal_fingerprint_distance(fp_1, fp_2, fp_method = fp_method))
    
    ### test get_scaffold
    # import sys
    #
    # path_list = sys.path
    # module_path = '/Users/guohan/Documents/Code/Tool'
    # if module_path not in sys.path:
    #     sys.path.append(module_path)
    #     print('Add module path')
    #
    # from process_SMILES import plot_2d_molecule_from_smiles
    #
    # smiles_1 = 'CCOC(=O)N1CCN(CC1)C1=CC=CC=C1NS(=O)(=O)C1=CC=C(C=C1)S(=O)(=O)N(C)C'
    # smiles_2 = 'COC(C)(C)C(=O)N1CCN(CC1)C1=C(NS(=O)(=O)C2=CC=C(C=C2)S(=O)(=O)N(C)C)C=CC=C1'
    #
    # plot_2d_molecule_from_smiles(smiles_1, output_file_without_ext = 'smiles_1_original', legend = 'smiles_1', atomNumber = None)
    # plot_2d_molecule_from_smiles(smiles_2, output_file_without_ext = 'smiles_2_original', legend = 'smiles_2', atomNumber = None)
    #
    # scaffold_smiles_1 = get_scaffold(smiles_1, include_chirality=False, generic = False)
    # scaffold_smiles_2 = get_scaffold(smiles_2, include_chirality=False, generic = False)
    #
    # plot_2d_molecule_from_smiles(scaffold_smiles_1, output_file_without_ext = 'smiles_1', legend = 'smiles_1', atomNumber = None)
    # plot_2d_molecule_from_smiles(scaffold_smiles_2, output_file_without_ext = 'smiles_2', legend = 'smiles_2', atomNumber = None)
    #
