#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 13:33:55 2021

@author: guohan

1. cluster SMILES by fingerprint
2. cluster SMILES by scaffold

"""

import sys, warnings

path_list = sys.path
module_path = '/Users/guohan/Documents/Code/Tool/utils'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from molecular_description import get_fingerprint, cal_fingerprint_distance, get_scaffold



########################
###### Clustering ######
########################

### Helper functions ###

def cluster_from_dist_matrix(dist_matrix, threshold):
    """
    cluster by DBSCAN, the similarity smaller than threshold is not clustered in one cluster,
    DBSCAN eps = 1.0 - threshold
    :para dist_matrix: np.ndarray matrix, similarity matrix
    :para threshold: float, similarity threshold, 1.0 - eps
    :return: list of ints, cluster labels
    """    
    # define DBSCAN object
    eps = 1.0 - threshold
    dbscan = DBSCAN(eps = eps, min_samples = 1, metric = "precomputed")
    
    # clustering
    labels = dbscan.fit_predict(dist_matrix)
    if False:
        print(labels)
    
    num_clusters = np.max(labels) + 1
    
    # treat each outlier as a cluster
    for i, label in enumerate(labels):
            
        if label == -1:
            labels[i] = num_clusters
            num_clusters += 1
        
    return labels


def labels_to_clusters(smiles_list, labels):
    """
    change a list of labels to clusters
    :para smiles_list: list of strs, list of SMILES
    :para labels: list of int, list of labels for SMILES
    :return: list of lists, clusters
    """
    assert len(smiles_list) == len(labels), 'Error: The smiles list and labels have different lengths'
    
    num_clusters = np.max(labels) + 1
    clusters = [[] for i in range(num_clusters)]
    
    for i, label in enumerate(labels):
        clusters[label].append(smiles_list[i])
    
    return clusters   


### Clustering by fingerprint ###

def get_fingerprint_dist_matrix(smiles_list, fp_method = 'ecfp4'):
    """
    compute distance matrix
    :para smiles_list: list of strs, list of SMILES
    :para fp_method: str, method to compute fingerprint (topology: topological fingerprint; 
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


def clustering_by_fingerprint(smiles_list, fp_method = 'ecfp4', cluster_threshold = 0.5, output_type = 'labels'):
    """
    cluster 'smiles_list' based on fingerprint with method 'fp_method',
    cluster by DBSCAN, the similarity smaller than threshold is not clustered in one cluster,
    DBSCAN eps = 1.0 - threshold
    :para smiles_list: list of strs, list of SMILES
    :para fp_method: str, method to compute fingerprint (topology: topological fingerprint; 
                                                         maccs: MACCS fingerprint;
                                                         atomPairs: atom pairs;
                                                         torsions: topological torsions;
                                                         ecfp4: Morgan ECFP4 (connectivity) fingerprint;
                                                         fcfp4: Morgan FCFP4 (feature-based) fingerprint)
                                                         mhfp6: MinHash fingerprint
    :para cluster_threshold: float, similarity threshold, 1.0 - eps
    :para output_type: str, output type, ('labels', 'clusters')
    :return: list of ints (labels) or list of lists (clusters)
    """
    smiles_list = np.array(smiles_list)
    fp_method = fp_method.lower()
    
    # calculate dist_matrix, while each elements is 1.0 - similarity
    dist_matrix = get_fingerprint_dist_matrix(smiles_list, fp_method)
    if False:
        df = pd.DataFrame(dist_matrix)
        df.to_csv('distance_matrix.csv')
    
    # cluster using DBSCAN
    labels = cluster_from_dist_matrix(dist_matrix, threshold = cluster_threshold)
    num_labels = np.max(labels) + 1
    
    print('Number of clusters:', num_labels)
    if True:
        for label in range(num_labels):
            print('Number elements in cluster {} is {}'.format(label, (labels == label).sum()))
    
    # return labels or clusters
    if output_type == 'labels':
        return labels
    elif output_type == 'clusters':
        clusters = labels_to_clusters(smiles_list, labels)
        return clusters
    else:
        warnings.warn('Error: Output type is not defined properly, return None')


### Clustering by scaffold ###

def clustering_by_scaffold(smiles_list, output_type = 'labels'):
    """
    cluster 'smiles_list' by scaffold
    :para smiles_list: list of strs, list of SMILES
    :para output_type: str, output type, ('labels', 'clusters')
    :return: list of ints (labels) or list of lists (clusters)
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
    if False:
        df = pd.DataFrame([i for i in range(num_scaffold)], columns = ['ID'])
        df['Scaffold'] = scaffolds
        df.to_csv('Scaffold.csv')
    
    # change to labels
    labels = [-1 for i in range(num_smiles)]
    
    for ID, scaffold in enumerate(scaffolds):
        for smiles_idx in clusters_dict[scaffold]:
            labels[smiles_idx] = ID
    
    labels = np.array(labels)

    print('Number of scaffolds:', num_scaffold)
    if True:
        for label in range(num_scaffold):
            print('Number elements in cluster {} is {}'.format(label, (labels == label).sum()))
    
    # return labels or clusters
    if output_type == 'labels':
        return labels
    elif output_type == 'clusters':
        clusters = labels_to_clusters(smiles_list, labels)
        return clusters
    else:
        warnings.warn('Error: Output type is not defined properly, return None')



if __name__ == '__main__':
    
    smiles_list = ['OC(CN1C=NC=N1)(CN1C=NC=N1)C1=C(F)C=C(F)C=C1', 'CCO', 'OCC', 'COO', 'c1ccccn1', 'c1ccco1',
                   'CCOC1=C(C=C(C=C1)S(=O)(=O)N(C)C)C2=NC(=O)C3=C(N2)C(=NN3C)C(C)(C)C',
                   'CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C']   
    
    print('*************************')
    print('Clustering by fingerprint')
    clusters = clustering_by_fingerprint(smiles_list, fp_method = 'ecfp4', cluster_threshold = 0.5, output_type = 'clusters')
    print(clusters)
    
    print('*************************')
    print('Clustering by scaffold')
    clusters = clustering_by_scaffold(smiles_list, output_type = 'clusters')
    print(clusters)
    print('*************************')
    