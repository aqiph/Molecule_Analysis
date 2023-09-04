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

from utils.clustering import fingerprint_based_clustering, mcs_based_clustering, scaffold_based_clustering, find_representative
from utils.tools import remove_unnamed_columns


def main(args):
    """
    group compounds in the args.input_file into clusters,
    based on fingerprint or scaffold
    """
    lprint = args.lprint
    input_file = os.path.abspath(args.input_file)

    ### clustering
    if args.lclustering:
        df = pd.read_csv(input_file)
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
        df['Cluster_label'] = labels

        df = df.reset_index(drop=True)
        df = remove_unnamed_columns(df)
        print('Number of rows in clustered file:', df.shape[0])
        output_file = f'{os.path.splitext(input_file)[0]}_clustered.csv'
        df.to_csv(output_file)
        input_file = output_file

    ### extracting representative compounds
    if args.lfinding_representative:
        df = pd.read_csv(input_file)
        COLUMNS = df.columns.tolist()
        df_representative = pd.DataFrame(columns = COLUMNS)

        label_list = df[args.label_column_name]
        for label in set(label_list):
            df_cluster = df[df[args.label_column_name] == label]
            df_cluster = df_cluster.reset_index(drop = True)
            smiles_list = df_cluster[args.smiles_column_name].values.tolist()
            indices = find_representative(smiles_list, args.representative_method, args.representative_counts)
            print(indices)
            df_new_representative = df_cluster.loc[indices]
            df_representative = pd.concat([df_representative, df_new_representative], ignore_index = True, sort = False)

        df_representative = df_representative.reset_index(drop = True)
        df_representative = remove_unnamed_columns(df_representative)
        print('Number of rows in representative file:', df_representative.shape[0])
        output_file = f'{os.path.splitext(input_file)[0]}_representative.csv'
        df_representative.to_csv(output_file)
    

def get_parser():
    """
    generate parser
    """
    argparser = argparse.ArgumentParser()
    
    argparser.add_argument('--input_file', default='examples.csv', type=str, help='Path of the input file.')
    argparser.add_argument('--lprint', default=False, action='store_true', help='Include this flag to enable printing distance matrix or scaffold.')
    argparser.add_argument('--smiles_column_name', default='Cleaned_SMILES', type=str, help='Name of the SMILES column.')

    ### clustering
    argparser.add_argument('--lclustering', default=False, action='store_true', help='Include this flag to enable clustring.')
    argparser.add_argument('--cluster_method', default='fingerprint', type=str, help='Method for compound clustering.')
    argparser.add_argument('--cluster_threshold', default=0.7, type=float, help='Similarity threshold for clustering.')
    argparser.add_argument('--fingerprint_method', default='ecfp4', type=str, help='Method for calculating fingerprint.')

    ### finding representative compounds
    argparser.add_argument('--lfinding_representative', default=False, action='store_true', help='Include this flag to enable finding representative compounds.')
    argparser.add_argument('--label_column_name', default='Cluster_label', type=str, help='Name of the label column.')
    argparser.add_argument('--representative_method', default='fingerprint', type=str, help='Method for finding representative compounds.')
    argparser.add_argument('--representative_counts', default=1, type=int, help='Number of representative compounds')


    args = argparser.parse_args()
    return args



if __name__ == '__main__':
    
    args = get_parser()
    main(args)
    
    
    