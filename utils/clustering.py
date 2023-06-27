#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 13:33:55 2021

@author: guohan

1. cluster compounds based on fingerprint
2. cluster compounds based on MCS
3. cluster compounds based on scaffold

"""

import sys

path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Molecule_Analysis/utils'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from molecular_description import get_fingerprint, cal_fingerprint_distance, cal_MCS, get_scaffold



########################
###### Clustering ######
########################

### Helper functions ###

def cluster_from_dist_matrix(dist_matrix, threshold):
    """
    cluster by DBSCAN, the similarity smaller than threshold is not clustered in one cluster,
    DBSCAN eps = 1.0 - threshold
    :param dist_matrix: np.ndarray matrix, similarity matrix
    :param threshold: float, similarity threshold, 1.0 - eps
    :return: list of ints, cluster labels
    """    
    # define DBSCAN object
    eps = 1.0 - threshold
    dbscan = DBSCAN(eps = eps, min_samples = 1, metric = "precomputed")
    
    # clustering
    labels = dbscan.fit_predict(dist_matrix)
    num_clusters = np.max(labels) + 1
    
    # treat each outlier as a cluster
    for i, label in enumerate(labels):
            
        if label == -1:
            labels[i] = num_clusters
            num_clusters += 1
        
    return labels


### Clustering based on fingerprint ###

def get_fingerprint_dist_matrix(smiles_list, fp_method = 'ecfp4'):
    """
    compute distance matrix based on fingerprint
    :param smiles_list: list of strs, list of SMILES
    :param fp_method: str, method to compute fingerprint (topology: topological fingerprint;
                                                         maccs: MACCS fingerprint;
                                                         atomPairs: atom pairs;
                                                         torsions: topological torsions;
                                                         ecfp4: Morgan ECFP4 (connectivity) fingerprint;
                                                         fcfp4: Morgan FCFP4 (feature-based) fingerprint)
                                                         mhfp6: MinHash fingerprint
    :return: np.ndarray matrix, similarity matrix
    """
    fp_method = fp_method.lower()
    
    # calculate fingerprints
    num_smiles = len(smiles_list)
    fp_list = [get_fingerprint(smiles, fp_method) for smiles in smiles_list]
    
    # calculate distance matrix as np.ndarray, distance = 1 - similarity
    dist_matrix = np.zeros((num_smiles, num_smiles), dtype = np.float16)
    
    for i in range(1, num_smiles):
        fp_1 = fp_list[i]
        for j in range(i):
            
            fp_2 = fp_list[j]
            dist = cal_fingerprint_distance(fp_1, fp_2, fp_method)
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    
    return np.array(dist_matrix)


def fingerprint_based_clustering(smiles_list, fp_method = 'ecfp4', cluster_threshold = 0.5, lprint =False):
    """
    cluster 'smiles_list' based on fingerprint with method 'fp_method',
    cluster by DBSCAN, the similarity smaller than threshold is not clustered in one cluster,
    DBSCAN eps = 1.0 - threshold
    :param smiles_list: list of strs, list of SMILES
    :param fp_method: str, method to compute fingerprint (topology: topological fingerprint;
                                                         maccs: MACCS fingerprint;
                                                         atomPairs: atom pairs;
                                                         torsions: topological torsions;
                                                         ecfp4: Morgan ECFP4 (connectivity) fingerprint;
                                                         fcfp4: Morgan FCFP4 (feature-based) fingerprint)
                                                         mhfp6: MinHash fingerprint
    :param cluster_threshold: float, similarity threshold, 1.0 - eps
    :param lprint: bool, whether to print distance matrix or not
    :return: list of ints (cluster labels)
    """
    smiles_list = np.array(smiles_list)
    fp_method = fp_method.lower()
    
    # calculate dist_matrix, while each element is 1.0 - similarity
    dist_matrix = get_fingerprint_dist_matrix(smiles_list, fp_method)
    if lprint:
        df = pd.DataFrame(dist_matrix)
        df.to_csv('distance_matrix.csv')
    
    # cluster using DBSCAN
    labels = cluster_from_dist_matrix(dist_matrix, threshold = cluster_threshold)
    num_labels = np.max(labels) + 1
    
    print('Number of clusters:', num_labels)
    for label in range(num_labels):
        print('Number of compounds in cluster {} is {}.'.format(label, (labels == label).sum()))
    
    # return labels or clusters
    return labels


### Clustering based on MCS ###

def get_mcs_dist_matrix(smiles_list):
    """
    compute distance matrix based on Maximum Common Substructure
    :param smiles_list: list of strs, list of SMILES
    """
    # calculate distance matrix as np.ndarray, distance = 1 - similarity
    num_smiles = len(smiles_list)
    dist_matrix = np.zeros((num_smiles, num_smiles), dtype=np.float16)

    for i in range(1, num_smiles):
        smiles_1 = smiles_list[i]
        for j in range(i):
            smiles_2 = smiles_list[j]
            dist = 1.0 - cal_MCS(smiles_1, smiles_2)
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    return np.array(dist_matrix)


def mcs_based_clustering(smiles_list, cluster_threshold = 0.5, lprint=False):
    """
    cluster 'smiles_list' based on MCS,
    cluster by DBSCAN, the similarity smaller than threshold is not clustered in one cluster,
    DBSCAN eps = 1.0 - threshold
    :param smiles_list: list of strs, list of SMILES
    :param cluster_threshold: float, similarity threshold, 1.0 - eps
    :param lprint: bool, whether to print distance matrix or not
    :return: list of ints (cluster labels)
    """
    smiles_list = np.array(smiles_list)

    # calculate dist_matrix, while each element is 1.0 - similarity
    dist_matrix = get_mcs_dist_matrix(smiles_list)
    if lprint:
        df = pd.DataFrame(dist_matrix)
        df.to_csv('distance_matrix.csv')

    # cluster using DBSCAN
    labels = cluster_from_dist_matrix(dist_matrix, threshold=cluster_threshold)
    num_labels = np.max(labels) + 1

    print('Number of clusters:', num_labels)
    for label in range(num_labels):
        print('Number of compounds in cluster {} is {}.'.format(label, (labels == label).sum()))

    # return labels or clusters
    return labels


### Clustering based on scaffold ###

def scaffold_based_clustering(smiles_list, lprint = False):
    """
    cluster 'smiles_list' based on scaffold
    :param smiles_list: list of strs, list of SMILES
    :param lprint: bool, whether to print scaffolds or not
    :return: list of ints (cluster labels)
    """
    clusters_dict = {} # {scaffold_smiles: [smiles_id]}
    num_smiles = len(smiles_list)
    
    # calculate scaffold cluster dictionary
    for i, smiles in enumerate(smiles_list):
        scaffold = get_scaffold(smiles, include_chirality = False, generic = False)
        
        if scaffold not in clusters_dict:
            clusters_dict[scaffold] = []
        clusters_dict[scaffold].append(i)
    
    num_scaffold = len(clusters_dict)
    scaffolds = clusters_dict.keys()
    
    # write scaffold table
    if lprint:
        df = pd.DataFrame([i for i in range(num_scaffold)], columns = ['ID'])
        df['Scaffold'] = scaffolds
        df.to_csv('Scaffold.csv')
    
    # assign cluster label (scaffold label) to each SMILES
    labels = [-1 for _ in range(num_smiles)]
    for ID, scaffold in enumerate(scaffolds):
        for smiles_idx in clusters_dict[scaffold]:
            labels[smiles_idx] = ID
    labels = np.array(labels)

    print('Number of scaffolds:', num_scaffold)
    for label in range(num_scaffold):
        print('Number of compounds in cluster {} is {}.'.format(label, (labels == label).sum()))
    
    # return labels or clusters
    return labels



if __name__ == '__main__':
    
    smiles_list = ['OC(CN1C=NC=N1)(CN1C=NC=N1)C1=C(F)C=C(F)C=C1', 'CCO', 'OCC', 'COO', 'c1ccccn1', 'c1ccco1',
                   'CCOC1=C(C=C(C=C1)S(=O)(=O)N(C)C)C2=NC(=O)C3=C(N2)C(=NN3C)C(C)(C)C',
                   'CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C']   
    
    print('*************************')
    print('Clustering based on fingerprint')
    clusters = fingerprint_based_clustering(smiles_list, fp_method = 'ecfp4', cluster_threshold = 0.5, lprint = False)
    print(clusters)

    print('*************************')
    print('Clustering based on MCS')
    clusters = mcs_based_clustering(smiles_list, cluster_threshold=0.5, lprint=False)
    print(clusters)
    
    print('*************************')
    print('Clustering based on scaffold')
    clusters = scaffold_based_clustering(smiles_list, lprint = False)
    print(clusters)
    print('*************************')
    