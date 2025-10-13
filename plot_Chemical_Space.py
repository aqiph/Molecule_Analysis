#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 14:00:00 2025

@author: guohan

"""


import os
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from utils.molecular_description import get_fingerprint
from collections import OrderedDict


def figure_plot(df, output_folder):
    sources = df['Source'].unique()
    data = OrderedDict()
    if 'Trainset' in sources:
        data['Trainset'] = df[df['Source'] == 'Trainset']
    if 'Trainset_positive' in sources:
        data['Trainset_negative'] = df[df['Source'] == 'Trainset_negative']
        data['Trainset_positive'] = df[df['Source'] == 'Trainset_positive']
    if 'Testset' in sources:
        # data['Testset'] = df[(df['Source'] == 'Testset')]
        data['Testset'] = df[(df['Source'] == 'Testset') & (df['Selected'] == 1)]

    plt.figure(figsize=(8, 8))
    color_list = ['r', 'b', 'orange', 'g']

    for i, (key, value) in enumerate(data.items()):
        if value is None:
            continue
        plt.scatter(value['tSNE1'], value['tSNE2'], color=color_list[i], label=key, s=50, alpha=0.8)

    plt.xlabel("t-SNE Dimension 1", fontsize=16)
    plt.ylabel("t-SNE Dimension 2", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('t-SNE Visualization of Chemical Space', fontsize=18, fontweight='bold')
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'Chemical_space.pdf'), format='pdf', bbox_inches='tight')


def plot_chemical_space_tSNE(input_file_library,
              input_file_query = None,
              smiles_column_name_library = 'SMILES',
              smiles_column_name_query = 'SMILES',
              label_column_name = None):
    """
    Plot t-SME for chemical space.
    :param input_file: str, path of the input file.
    :param smiles_column_name: str, name of the SMILES column.
    :param label_column_name: str, name of the label column.
    """
    # get library data (training data)
    df_library = pd.read_csv(input_file_library)
    df_library.rename(columns = {smiles_column_name_library : 'SMILES'}, inplace = True)
    if label_column_name is None:
        df_library['Source'] = 'Trainset'
    else:
        df_library['Source'] = df_library[label_column_name].apply(lambda l: 'Trainset_positive' if l >= 1 else 'Trainset_negative')

    # get query data (test data)
    if input_file_query is None:
        df_query = None
    else:
        df_query = pd.read_csv(input_file_query)
        df_query.rename(columns = {smiles_column_name_query: 'SMILES'}, inplace = True)
        df_query['Source'] = 'Testset'

    # combine
    if df_query is not None:
        df_combined = pd.concat([df_library, df_query], ignore_index=True)
    else:
        df_combined = df_library

    # calculate fingerprint
    df_combined['Fingerprint'] = df_combined['SMILES'].apply(get_fingerprint)
    df_combined['Fingerprint'] = df_combined['Fingerprint'].apply(lambda fp: np.array([int(bit) for bit in fp.ToBitString()], dtype=int))
    df_combined = df_combined.dropna(subset=['Fingerprint'])

    # Apply t-SNE to the combined fingerprints
    print("Running t-SNE dimensionality reduction...")
    fps = np.stack(df_combined['Fingerprint'].values)
    tsne = TSNE(n_components=2,
                perplexity=30,
                n_iter=2000,
                random_state=42,
                verbose=1
                )
    tsne_results = tsne.fit_transform(fps)
    df_combined['tSNE1'] = tsne_results[:, 0]
    df_combined['tSNE2'] = tsne_results[:, 1]

    # plot
    figure_plot(df_combined, output_folder=os.path.split(input_file_library)[0])



if __name__ == "__main__":

    input_file_library = 'plot_Chemical_Space/tests/test_plot_Chemical_Space_train.csv'
    input_file_query = 'plot_Chemical_Space/tests/test_plot_Chemical_Space_test.csv'
    plot_chemical_space_tSNE(input_file_library,
                             input_file_query=input_file_query,
                             smiles_column_name_library='SMILES',
                             smiles_column_name_query='SMILES',
                             label_column_name=None)

