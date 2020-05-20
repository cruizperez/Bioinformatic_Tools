#!/usr/bin/env python

"""
########################################################################
# Author:   Carlos Ruiz
# Email:    cruizperez3@gatech.edu
# Version:  1.0
# Date:     September 08 2019

# Description: This script transforms a tab-separated table into a matrix.
# You can select which columns will become rows, columns and values.
########################################################################
"""
################################################################################
"""---0.0 Import Modules---"""
import pandas as pd

################################################################################
"""---1.0 Define Functions---"""

def table_converter(table_file, outfile, first_column, second_column, value_col):
    row_ids = []
    col_ids = []
    
    with open(table_file, 'r') as input_file:
        for line in input_file:
            line = line.strip().split()
            row_ids.append(line[first_column-1])
            col_ids.append(line[second_column-1])

    final_matrix = pd.DataFrame(index=row_ids, columns=col_ids)

    with open(table_file, 'r') as input_file:
        for line in input_file:
            line = line.strip().split()
            final_matrix.loc[first_column-1, second_column-1] = value_col-1

    final_matrix.to_csv(outfile, sep="\t", header=True, index=True)

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script transforms a tab-separated table into a matrix.\n'''
            '''You can select which columns will become rows, columns and values.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [input table] -o [output matrix]\n'''
            '''Global mandatory parameters: -i [input table] -o [output matrix]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_table', action='store', required=True,
                        help='Input tab-delimited table to parse.')
    parser.add_argument('-o', '--output', dest='output_matrix', action='store', required=True,
                        help='Output matrix.')
    parser.add_argument('--rows', dest='rows', action='store', required=False, type=int, default=1,
                        help='Column in the input file that will make the rows of the matrix. By default the first (1).')
    parser.add_argument('--cols', dest='cols', action='store', required=False, type=int, default=2,
                        help='Column in the input file that will make the columns of the matrix. By default the second (2).')
    parser.add_argument('--values', dest='values', action='store', required=False, type=int, default=3,
                        help='Column in the input file with values to fill the matrix. By default the third (3).')
    args = parser.parse_args()

    input_table = args.input_table
    output_matrix = args.output_matrix
    rows = args.rows
    cols = args.cols
    values = args.values

    table_converter(input_table, output_matrix, rows, cols, values)

if __name__ == "__main__":
    main()