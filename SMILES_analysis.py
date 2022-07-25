#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 17:00:16 2021

@author: guohan

1. Group molecules into clusters based on fingerprint or scaffold
2. Determine whether a molecule is a neighbor of centroid molecule based on fingerprint

"""

import os
import pandas as pd
import argparse
from utils.clustering import fingerprint_clustering, scaffold_clustering
from utils.nearest_neighbor import fingerprint_nearest_neighbor


def cluster_analysis(args):
    """
    group compounds in the args.input_file into clusters,
    based on fingerprint or scaffold
    """
    # read files
    df = pd.read_csv(args.input_file, index_col = 0)
    
    # output name
    output_file_without_ext = os.path.splitext(os.path.abspath(args.input_file))[0]
    output_file = output_file_without_ext + '_cluster.csv'
    
    # clustering
    smiles_list = df['Cleaned_SMILES'].values.tolist()
    
    if args.cluster_method == 'fingerprint':
        fp_method = args.fingerprint_method
        cluster_threshold = args.cluster_threshold
        
        labels = fingerprint_clustering(smiles_list, fp_method, cluster_threshold, output_type = 'labels')
        
    elif args.cluster_method == 'scaffold':
        labels = scaffold_clustering(smiles_list, output_type = 'labels')
    
    # add cluster labels
    df['Cluster_label'] = labels
    
    # write to file
    df = df.reset_index(drop = True)
    
    print('Number of rows:', df.shape[0])
    df.to_csv(output_file)


def nearest_neighbor_analysis(args):
    """
    determine whether a molecule in df_candidate is a neighbor of molecules in df_centroid,
    based on fingerprint with method 'fp_method'.
    if the similarity is larger than neighbor_threshold,
    the candidate molecule is considered as a neighbor of centroid molecule
    """
    # read files
    df_centroid = pd.read_csv(args.input_file_centroid, index_col = 0)
    df_candidate = pd.read_csv(args.input_file, index_col = 0)
    
    # output name
    output_file_without_ext = os.path.splitext(os.path.abspath(args.input_file))[0]
    output_file = output_file_without_ext + '_neighbor.csv'
    
    centroid_smiles_list = df_centroid['Cleaned_SMILES'].values.tolist()
    candidate_smiles_list = df_candidate['Cleaned_SMILES'].values.tolist()
    
    # determine whether a candidate molecule is a neighbor of centroid molecules
    fp_method = args.fingerprint_method
    neighbor_threshold = float(args.neighbor_threshold)
    
    labels = fingerprint_nearest_neighbor(centroid_smiles_list, candidate_smiles_list, fp_method, neighbor_threshold)
    labels = [int(label) for label in labels]
    
    # add cluster labels
    df_candidate['Cluster_label'] = labels
    
    # write to file
    df_candidate = df_candidate.reset_index(drop = True)
    
    print('Number of rows:', df_candidate.shape[0])
    df_candidate.to_csv(output_file)
    

def get_parser():
    """
    generate parser
    """
    argparser = argparse.ArgumentParser()
    
    argparser.add_argument('-input_file', default = 'SMILES_analysis/tests/examples.csv', help = 'Input file name')
    argparser.add_argument('-input_file_centroid', default = 'SMILES_analysis/tests/examples_centroid.csv', help = 'Input file name for centroid data points')
    
    argparser.add_argument('-cluster_method', default = 'fingerprint', help = 'The method to cluster compound', type = str)    
    argparser.add_argument('-fingerprint_method', default = 'ecfp4', help = 'The method to calculate fingerprint', type = str)
    argparser.add_argument('-cluster_threshold', default = 0.5, help = 'The similarity threshold for clustering', type = float)
    argparser.add_argument('-neighbor_threshold', default = 0.5, help = 'The similarity thresohld for defining nearest neighbor data points')

    args = argparser.parse_args()
    return args



if __name__ == '__main__':
    
    args = get_parser()
    
    cluster_analysis(args)
    nearest_neighbor_analysis(args)
    
    
    