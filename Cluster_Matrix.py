#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 24 March 2019

# Description: This script clusters a matrix, plots a heatmap and outputs the original table reordered.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import argparse, sys
# Change default font used by matplotlib so it can be modified
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

################################################################################
"""---2.0 Define Functions---"""

def Matrix_Clustering(Input_file, Output_Prefix):
    Matrix = pd.read_csv(Input_file, sep="\t", index_col=0)
    Heatmap = sns.clustermap(Matrix, metric='euclidean', method='average', figsize= (20,20))
    Heatmap.savefig(Output_Prefix + ".pdf")

    # Extract reordered indices and colums
    Matrix_index = Heatmap.dendrogram_row.reordered_ind
    Matrix_cols = Heatmap.dendrogram_col.reordered_ind
    Index = Matrix.index[Matrix_index]
    Cols = Matrix.index[Matrix_cols]
    Reordered_Matrix = pd.DataFrame(Matrix, Index, columns=Cols)

    return Reordered_Matrix

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Clusters a matrix, plots a heatmap of the clustered matrix and outputs the original table reordered'''
                                    'Global mandatory parameters: [Matrix]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--inputMatrix", dest='Input_Matrix', required=True, help="Input Matrix file")
    parser.add_argument('-p', '--prefix', dest='Prefix', action='store', help='Output prefix, by default "Cluster_Matrix"')
    args = parser.parse_args()

    Input_Matrix = args.Input_Matrix
    Prefix = args.Prefix

    if Prefix == None:
        Prefix = "Cluster_Matrix"

    # Run cluster Function
    Clustered_Matrix = Matrix_Clustering(Input_Matrix, Prefix)
    Clustered_Matrix.to_csv(Prefix + ".tab", sep="\t")

if __name__ == "__main__":
    main()
