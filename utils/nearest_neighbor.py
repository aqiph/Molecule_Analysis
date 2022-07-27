#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 15:58:36 2022

@author: guohan

1. find the nearest neighbor of centroid SMILES from candidate SMILES

"""

import sys

path_list = sys.path
module_path = '/Users/guohan/Documents/Code/Tool/utils'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

from molecular_description import get_fingerprint, cal_fingerprint_distance



##############################
###### Nearest Neighbor ######
##############################

### Helper functions ###


### Find neighbors of centroid data points ###

def nearest_neighbor_by_fingerprint(centroid_smiles_list, candidate_smiles_list, fp_method = 'ecfp4', neighbor_threshold = 0.5):
    """
    determine whether a molecule in df_candidate is a neighbor of molecules in df_centroid,
    based on fingerprint with method 'fp_method'.
    if the similarity is larger than neighbor_threshold,
    the candidate molecule is considered as a neighbor of centroid molecule
    :para centroid_smiles_list: list of SMILES, the centroid molecules
    :para candidate_smiles_list: list of SMILES, the candidate molecules
    :para fp_method: str, method to compute fingerprint (topology: topological fingerprint; 
                                                         maccs: MACCS fingerprint;
                                                         atomPairs: atom pairs;
                                                         torsions: topological torsions;
                                                         ecfp4: Morgan ECFP4 (connectivity) fingerprint;
                                                         fcfp4: Morgan FCFP4 (feature-based) fingerprint)
                                                         mhfp6: MinHash fingerprint
    :para neighbor_threshold: float, similarity threshold
    """
    fp_method = fp_method.lower()
    
    # calculate fingerprints    
    centroid_fp_list = [get_fingerprint(smiles, fp_method) for smiles in centroid_smiles_list]
    candidate_fp_list = [get_fingerprint(smiles, fp_method) for smiles in candidate_smiles_list]
    
    num_centroid_smiles = len(centroid_smiles_list)
    num_candidate_smiles = len(candidate_smiles_list)
    
    labels = [False for _ in range(num_candidate_smiles)]
    
    # determine whether a candidate molecule is a neighbor of a centroid molecule
    for i in range(num_candidate_smiles):
        candidate_fp = candidate_fp_list[i]
        
        j = 0
        while (not labels[i]) and (j < num_centroid_smiles):
            centroid_fp = centroid_fp_list[j]
            
            dist = cal_fingerprint_distance(candidate_fp, centroid_fp, fp_method)
            sim = 1.0 - dist
            if sim >= neighbor_threshold:
                print('centroid {} candidate {}, sim {}'.format(j, i, str(sim)))
                labels[i] = True
            
            j += 1
    
    return labels



if __name__ == '__main__':
    
    centroid_smiles_list = ['OC(CN1C=NC=N1)(CN1C=NC=N1)C1=C(F)C=C(F)C=C1', 'CCO']
    candidate_smiles_list = ['OC(CN1C=NC=N1)(CN1C=NC=N1)C1=C(F)C=C(F)C=C1', 'CCO', 'OCC', 'COO', 'c1ccccn1', 'c1ccco1',
                   'CCOC1=C(C=C(C=C1)S(=O)(=O)N(C)C)C2=NC(=O)C3=C(N2)C(=NN3C)C(C)(C)C',
                   'CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C']
    
    labels = nearest_neighbor_by_fingerprint(centroid_smiles_list, candidate_smiles_list, fp_method = 'ecfp4', neighbor_threshold = 0.5)
    print(labels)
    



    
    
    


