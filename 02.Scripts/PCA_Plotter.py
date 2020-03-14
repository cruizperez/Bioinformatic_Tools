#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      0.9
# Date:         March 07, 2020

# Description: Normalizes and plots the PCA ordination of a given table.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
from numpy.random import randint

################################################################################
"""---1.0 Define Functions---"""

def table_normalizer_pca(input_table, headers=True, index=True, features=None, target=None):
    if index == False:
        index = None
    else:
        index = 0
    if headers == True:
        headers = 'infer'
    original_table = pd.read_csv(input_table, sep="\t", header=headers, index_col=index)
    print(original_table)
    # Extract target and features if necessary
    if target is not None:
        target_table = original_table.loc[:, target].values
    if features is not None:
        feature_table = original_table.loc[:, features].values
        feature_table_standard = StandardScaler().fit_transform(feature_table)
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(feature_table_standard)
        components_table = pd.DataFrame(principalComponents, columns=['PC1', 'PC2'])
        if features is not None and target is not None:
            final_table = pd.concat([components_table, target_table], axis=1)
    else:
        table_standard = StandardScaler().fit_transform(original_table)
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(table_standard)
        final_table = pd.DataFrame(principalComponents, columns=['PC1', 'PC2'])
    
    return final_table

def pca_plotter(input_table, prefix, target=None):
    matplotlib.rcParams['pdf.fonttype'] = 42
    if target is not None:
        targets = input_table.loc[:, target].values
        target_values = list(targets.unique)
        colors = []
        for _ in range(len(target_values)):
            colors.append(tuple(randint(256, size=3) + (255,))) 
        # Plot PCA
        Figure, Axis =  plt.subplots(1,1, figsize=(10,10), dpi=300, constrained_layout=True)
        Axis.set_xlabel('Principal Component 1', fontsize = 15)
        Axis.set_ylabel('Principal Component 2', fontsize = 15)
        Axis.set_title('2-Component PCA', fontsize = 20)
        for target, color in zip(target_values,colors):
            indicesToKeep = input_table[target] == target
            Axis.scatter(input_table.loc[indicesToKeep, 'PC1'], 
                        input_table.loc[indicesToKeep, 'PC2'],
                        c = color, s = 50)
        Axis.legend(targets)
        Figure.savefig(prefix + ".pdf", transparent=True)
    else:
        Figure, Axis =  plt.subplots(1,1, figsize=(10,10), dpi=300, constrained_layout=True)
        Axis.set_xlabel('Principal Component 1', fontsize = 15)
        Axis.set_ylabel('Principal Component 2', fontsize = 15)
        Axis.set_title('2-Component PCA', fontsize = 20)
        Axis.scatter(input_table.loc[:, 'PC1'], input_table.loc[:, 'PC2'], s = 50)
        Figure.savefig(prefix + ".pdf", transparent=True)

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''Normalizes and plots the PCA ordination of a given table.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input table] -p [Output Prefix]\n'''
            '''Global mandatory parameters: -i [Input table] -p [Output Prefix]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_table', action='store', required=True,
                        help='Input table')
    parser.add_argument('-p', '--prefix', dest='prefix', action='store', required=True,
                        help='Output prefix for normalized table and plot')
    parser.add_argument('--headers', dest='headers', action='store_false', required=False,
                        help='Input table has headers. By default True')
    parser.add_argument('--index', dest='index', action='store_false', required=False,
                        help='Input table has index. By default True')
    parser.add_argument('--features', dest='features', action='store', required=False, nargs='+',
                        help='List of colums to use as features. By default all.')
    parser.add_argument('--target', dest='target', action='store', required=False,
                        help='Column to use as target for grouping. By default None.')
    args = parser.parse_args()

    input_table = args.input_table
    prefix = args.prefix
    headers = args.headers
    index = args.index
    features = args.features
    target = args.target

    # ----------------------------
    table_for_pca = table_normalizer_pca(input_table, headers, index, features, target)
    pca_plotter(table_for_pca, prefix, target)
    # ----------------------------

if __name__ == "__main__":
    main()

    