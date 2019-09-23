#!/usr/bin/env python

"""
########################################################################
# Author:   Carlos Ruiz
# Email:    cruizperez3@gatech.edu
# Version:  1.0
# Date:     September 08 2019

# Description: This script merges tables based on a common ID.
# You can pass several tables and specify the column of the ID and columns to merge.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Table_Merger(Table_Files, Output, ID_Column=1, Item_Column=2, Header=False, Col_Names = None):
    import pandas as pd
    Table_List = []
    for Index, Table in enumerate(Table_Files):
        ID = ID_Column[Index].split(",")
        ID[:] = [int(x) - 1 for x in ID]
        Columns = Item_Column[Index].split(",")
        Columns[:] = [int(x) - 1 for x in Columns]
        Total_Cols = ID + Columns
        if Header == True:
            New_Table = pd.read_csv(Table, sep = "\t", header=0, usecols = Total_Cols, index_col=ID)
        else:
            New_Table = pd.read_csv(Table, sep = "\t", header=None, usecols = Total_Cols, index_col=ID)
        Table_List.append(New_Table)

    Final_Table = pd.concat(Table_List, axis=1, sort=False)
    if Col_Names != None:
        Final_Table.index.name = Col_Names[0]
        Final_Table.columns = Col_Names[1:]
    Final_Table.to_csv(Output, sep="\t")






################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script builds a sqlite database from a tab separated table.\n'''
            '''By default it assumes the first line of the input table has headers, if not\n'''
            '''you must provide a list of headers as -h header1,header2,header3...\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input Table] -d [Database Name] -t [Table Name] -n [Index Name] --header [Header List]\n'''
            '''Global mandatory parameters: -i [Input Table] -d [Database Name] -t [Table Name] -n [Index Name]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-t', '--tables', dest='Input_Tables', action='store', nargs='+', required=True, help='Input tab-delimited table to parse, by default assumes headers are present')
    parser.add_argument('-o', '--output', dest='Output_Table', action='store', required=True, help='Input tab-delimited table to parse, by default assumes headers are present')
    parser.add_argument('--id', dest='ID_Locations', action='store', nargs='+', required=False, help='Output tab-delimited table to store annotations')
    parser.add_argument('--cols', dest='Columns', action='store', nargs='+', required=True, help='SQL database where the annotations are stored')
    parser.add_argument('--header', dest='Header', action='store_true', required=False, help='Column with gene IDs in the input file, by default None, i.e. assumes only IDs of hits are given.')
    parser.add_argument('--name_col', dest='Col_Names', action='store', nargs='+', required=False, help='Column with IDs of database hits in the input file')
    args = parser.parse_args()

    Input_Tables = args.Input_Tables
    Output_Table = args.Output_Table
    ID_Locations = args.ID_Locations
    Columns = args.Columns
    Header = args.Header
    Col_Names = args.Col_Names

    Table_Merger(Input_Tables, Output_Table, ID_Locations, Columns, Header, Col_Names)

if __name__ == "__main__":
    main()