#!/usr/bin/env python

# ============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Jan 27th, 2021

# Description: Creates a heatmap from an input table. It 
# and calculates the completeness percentage of each module present.
# ============================================================================


# ============================================================================
# Import required modules
import argparse
import seaborn as sns
import pandas as pd
import matplotlib
from sys import argv
from typing import Any
matplotlib.rcParams['pdf.fonttype'] = 42
# ============================================================================


# ============================================================================
# Define functions
# Read table
def read_table(input_table: str, headers: Any,
               rows: Any, transpose: bool) -> pd.DataFrame:
    input_dataframe = pd.read_csv(input_table, sep="\t", 
                                  header=headers, index_col=rows)
    if transpose == True:
        input_dataframe = input_dataframe.T

    return input_dataframe

# Plot heatmap
def heatmap_plot(input_table: pd.DataFrame, output: str, bar_label: str, 
                 xlabel: str, ylabel: str, cluster_rows: bool,
                 cluster_cols: bool) -> None:
    sns.set(font_scale=1.6)
    heatmap = sns.clustermap(input_table, xticklabels=True, yticklabels=True, 
                         figsize=(10,10), row_cluster=cluster_rows, 
                         col_cluster=cluster_cols, 
                         cbar_pos=(0.01, 1, 0.2, 0.04),
                         cbar_kws= {'orientation': 'horizontal', 
                                    'label': f"{bar_label}"},
                         dendrogram_ratio=0.1)
    heatmap.ax_heatmap.set_xlabel(f"{xlabel}", fontsize=20)
    heatmap.ax_heatmap.set_ylabel(f"{ylabel}", fontsize=20)
    heatmap.savefig(output)
# ============================================================================


# ============================================================================
# Define main function
def main():
    # Setup parser for arguments
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description='''Reads a tab-separated table and plots a heatmap\n'''
            '''Usage: ''' + argv[0] + ''' -i [Input Table] -o [Output]\n'''
            '''Global mandatory parameters: -i [Input Table] -o [Output]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    mandatory_arguments = parser.add_argument_group('Required')
    mandatory_arguments.add_argument('-i', '--input_table', dest='input_table',
                        action='store', required=True,
                        help='Input tab-separated table.')
    mandatory_arguments.add_argument('-o', '--output', dest='output', 
                        action='store', required=True,
                        help='Output file, can have .png or .pdf extension.')
    
    optional_arguments = parser.add_argument_group('optional')
    optional_arguments.add_argument('--headers', dest='headers', 
                        action='store', required=False, default=None, type=int,
                        help='Row to use as table header. Default no header.')
    optional_arguments.add_argument('--rows', dest='rows', 
                        action='store', required=False, default=None, type=int,
                        help='Col to use as row names. Default no row names.')
    optional_arguments.add_argument('--transpose', dest='transpose', 
                        action='store_true', required=False,
                        help='Transpose input table.')
    optional_arguments.add_argument('--bar_label', dest='bar_label', 
                        action='store', required=False, default=None, type=str,
                        help='Label of heatmap color bar. Default no label.')
    optional_arguments.add_argument('--xlabel', dest='xlabel', 
                        action='store', required=False, default=None, type=str,
                        help='Label of X axis. Default, no label.')
    optional_arguments.add_argument('--ylabel', dest='ylabel', 
                        action='store', required=False, default=None, type=str,
                        help='Label of Y axis. Default, no label.')
    optional_arguments.add_argument('--cluster_rows', dest='cluster_rows', 
                        action='store_true', required=False, 
                        help='Cluster heatmap rows. Default, no clustering.')
    optional_arguments.add_argument('--cluster_cols', dest='cluster_cols', 
                        action='store_true', required=False, 
                        help='Cluster heatmap cols. Default, no clustering.')
    # If no arguments provided
    if len(argv) == 1:
        parser.print_help()
        exit(1)
    
    # Parse variables
    args = parser.parse_args()
    input_table = args.input_table
    output = args.output
    headers = args.headers
    if headers is not None:
        headers -= 1
    rows = args.rows
    if rows is not None:
        rows -= 1
    transpose = args.transpose
    bar_label = args.bar_label
    if bar_label is None:
        bar_label = ""
    xlabel = args.xlabel
    ylabel = args.ylabel
    cluster_rows = args.cluster_rows
    cluster_cols = args.cluster_cols


    # Read input table
    data_table = read_table(input_table, headers, rows, transpose)

    # Plot heatmap
    heatmap_plot(data_table, output, bar_label, xlabel, ylabel, cluster_rows,
                 cluster_cols)
# ============================================================================


# ============================================================================
# Run main if executed from cli
if __name__ == '__main__':
    main()
# ============================================================================