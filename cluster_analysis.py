#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 17:00:16 2021

@author: guohan

1. Group molecules into clusters based on fingerprint or scaffold
2. Determine whether a molecule is a neighbor of centroid molecule based on fingerprint

"""

import os
import argparse
import pandas as pd

from utils.clustering import fingerprint_based_clustering, mcs_based_clustering, scaffold_based_clustering
from utils.tools import remove_unnamed_columns


def main(args):
    """
    group compounds in the args.input_file into clusters,
    based on fingerprint or scaffold
    """
    # read files
    input_file = os.path.abspath(args.input_file)
    df = pd.read_csv(input_file)
    
    # output file
    output_file_without_ext = os.path.splitext(input_file)[0]
    output_file = output_file_without_ext + '_clustered.csv'
    lprint = args.lprint
    
    # clustering
    smiles_list = df[args.smiles_column_name].values.tolist()
    
    if args.cluster_method == 'fingerprint':
        fp_method = args.fingerprint_method
        cluster_threshold = args.cluster_threshold
        labels = fingerprint_based_clustering(smiles_list, fp_method, cluster_threshold, lprint)

    elif args.cluster_method == 'mcs':
        cluster_threshold = args.cluster_threshold
        labels = mcs_based_clustering(smiles_list, cluster_threshold, lprint)
        
    elif args.cluster_method == 'scaffold':
        labels = scaffold_based_clustering(smiles_list, lprint)

    else:
        raise Exception('Error: Invalid cluster method.')
    
    # add cluster labels
    df['Cluster_label'] = labels
    
    # write to file
    df = df.reset_index(drop = True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)
    

def get_parser():
    """
    generate parser
    """
    argparser = argparse.ArgumentParser()
    
    argparser.add_argument('--input_file', default='examples.csv', type=str, help='Path of the input file.')
    argparser.add_argument('--smiles_column_name', default='Cleaned_SMILES', type=str, help='Name of the SMILES column.')
    argparser.add_argument('--lprint', default=False, action='store_true', help='Include this flag to enable printing distance matrix or scaffold.')

    argparser.add_argument('--cluster_method', default='fingerprint', type=str, help='The method to cluster compound.')
    argparser.add_argument('--fingerprint_method', default='ecfp4', type=str, help='The method to calculate fingerprint.')
    argparser.add_argument('--cluster_threshold', default=0.7, type=float, help='The similarity threshold for clustering.')

    args = argparser.parse_args()
    return args



if __name__ == '__main__':
    
    args = get_parser()
    main(args)
    
    
    