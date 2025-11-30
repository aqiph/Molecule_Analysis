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


def figure_plot(df, output_folder, plot_order=None, color_order=None):
    if plot_order is not None:
        sources = plot_order.copy()
    else:
        sources = df['Source'].unique()

    data = OrderedDict()
    for source in sources:
        data[source] = df[df['Source'] == source]

    plt.figure(figsize=(8, 8))
    if color_order is not None:
        color_list = color_order.copy()
        color_list += ['b', 'orange', 'r', 'g', 'purple', 'grey'][len(color_order):6]
    else:
        color_list = ['b', 'orange', 'r', 'g', 'purple', 'grey']

    for i, (key, value) in enumerate(data.items()):
        if value is None:
            continue
        plt.scatter(value['tSNE1'], value['tSNE2'], color=color_list[i%6], label=key, s=50, alpha=0.8)

    plt.xlabel("t-SNE Dimension 1", fontsize=16)
    plt.ylabel("t-SNE Dimension 2", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.xlim(-150, 150)
    # plt.ylim(-150, 150)
    plt.title('t-SNE Visualization of Chemical Space', fontsize=18, fontweight='bold')
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'Chemical_space.pdf'), format='pdf', bbox_inches='tight')


def plot_chemical_space_tSNE(
        input_file_library,
        input_file_query = None,
        smiles_column_name_library = 'SMILES',
        smiles_column_name_query = 'SMILES',
        label_column_name_library = None,
        label_column_name_query = None,
        plot_order = None,
        color_order = None):
    """
    Plot t-SME for chemical space.
    :param input_file_library: str, path of the input library file.
    :param input_file_query: str, path of the input query file.
    :param smiles_column_name_library: str, name of the SMILES column in the input library file.
    :param smiles_column_name_query: str, name of the SMILES column in the input query file.
    :param label_column_name_library: str, name of the label column in the input library file.
    :param label_column_name_query: str, name of the label column in the input query file.
    :param plot_order: list, order of the dots in the chemical space plot.
    :param color_order: list, order of the colors in the chemical space plot.
    """
    # get library data (training data)
    df_library = pd.read_csv(input_file_library)
    df_library.rename(columns = {smiles_column_name_library : 'SMILES'}, inplace = True)
    if label_column_name_library is None:
        df_library['Source'] = 'Trainset'
    else:
        df_library['Source'] = df_library[label_column_name_library]

    # get query data (test data)
    if input_file_query is None:
        df_query = None
    else:
        df_query = pd.read_csv(input_file_query)
        df_query.rename(columns = {smiles_column_name_query: 'SMILES'}, inplace = True)
        if label_column_name_query is None:
            df_query['Source'] = 'Testset'
        else:
            df_query['Source'] = df_query[label_column_name_query]

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
    figure_plot(df_combined, output_folder=os.path.split(input_file_library)[0], plot_order=plot_order, color_order=color_order)



if __name__ == "__main__":

    input_file_library = 'plot_Chemical_Space/tests/test_plot_Chemical_Space_train.csv'
    input_file_query = 'plot_Chemical_Space/tests/test_plot_Chemical_Space_test.csv'
    plot_chemical_space_tSNE(input_file_library,
                             input_file_query=input_file_query,
                             smiles_column_name_library='SMILES',
                             smiles_column_name_query='SMILES',
                             label_column_name_library=None,
                             label_column_name_query=None,
                             plot_order=None,
                             color_order=None)


