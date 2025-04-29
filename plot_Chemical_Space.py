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


def plot_tSNE(input_file, smiles_column_name = 'SMILES', label_column_name = 'Label'):
    """
    Plot t-SME for chemical space.
    :param input_file: str, path of the input file.
    :param smiles_column_name: str, name of the SMILES column.
    :param label_column_name: str, name of the label column.
    """
    # get training data
    df = pd.read_csv(input_file)
    df.sort_values(by=label_column_name, ascending=True, inplace=True)
    SMILESs = df[smiles_column_name].tolist()
    value_counts = df[label_column_name].value_counts()
    unique_labels = value_counts.index.tolist()
    unique_label_counts = value_counts.values.tolist()
    num_clusters = len(unique_label_counts)

    # get fingerprint
    fps = []
    for smiles in SMILESs:
        fp = get_fingerprint(smiles)
        fp = np.array([int(bit) for bit in fp.ToBitString()], dtype=int)
        fps.append(fp)
    fps = np.array(fps)

    # Apply t-SNE to the combined fingerprints
    tsne = TSNE(n_components=num_clusters, random_state=669, perplexity=50)
    tsne_results = tsne.fit_transform(fps)

    tsne_lists = []
    start_index = 0
    for num in unique_label_counts:
        end_index = start_index + num
        tsne_list = tsne_results[start_index:end_index]
        tsne_lists.append(tsne_list)
        start_index = end_index

    # plot
    plt.figure(figsize=(8, 8))
    color_list = ['r', 'b', 'orange', 'g']
    for i, tsne_list in enumerate(tsne_lists):
        plt.scatter(tsne_list[:, 0], tsne_list[:, 1], color = color_list[i%len(color_list)], label=f'{unique_labels[i]}', alpha=0.7)

    plt.xlabel("t-SNE Component 1", fontsize=14)
    plt.ylabel("t-SNE Component 2", fontsize=14)
    plt.xticks([])
    plt.yticks([])
    plt.title('Chemical Space', fontsize=16, fontweight='bold')
    plt.legend(fontsize=14)
    plt.tight_layout()
    output_folder = os.path.split(input_file)[0]
    plt.savefig(os.path.join(output_folder, 'Chemical_space.pdf'), format='pdf')



if __name__ == "__main__":

    input_file = 'plot_Chemical_Space/tests/test_chemical_space.csv'
    plot_tSNE(input_file)

