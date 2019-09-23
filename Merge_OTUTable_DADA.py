#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 August 19 2019

# Description: This script merges OTU tables from DADA2 into a single one.
# It takes a representative fasta sequence and multiple OTU tables to be merged.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def OTUTable_Merger(Input_Fasta, OTU_Table_List, Transpose = False):
    import pandas as pd
    Representative_Sequences = []
    Dataframe_List = []
    Input_Tables = OTU_Table_List
    
    # Parse FastA file and append IDs to Representative Sequences list.
    with open(Input_Fasta) as Fasta_Input:
        for line in Fasta_Input:
            if ">" in line:
                line = line.strip()
                line = line.replace(">", "")
                Representative_Sequences.append(line)

    # Sort the IDs of representative sequences by length
    Representative_Sequences.sort(key=len, reverse = True)

    # Open all tables, convert them to dataframes and add them to the list
    for table in Input_Tables:
        Dataframe_List.append(pd.read_csv(table, header = 0, 
                  sep = "\t", index_col = 0))

    # Concatenate OTU tables
    Merged_Table = pd.concat(Dataframe_List, axis=1, sort=True)
    Merged_Table = Merged_Table.fillna(0)
    Merged_Table_Final = Merged_Table.copy(deep=True)

    # Add IDs not found in the table to a temporal list
    # to reduce the number of columns to search for.
    Temporal_IDs = []

    for Column in Merged_Table:
        if Column in Representative_Sequences:
            continue
        else:
            Temporal_IDs.append(Column)

    # Search for substrings in the representative sequences.
    # Rename substrings by the long names.
    for Column in Temporal_IDs:
        for RepSeq in Representative_Sequences:
            if Column in RepSeq:
                Merged_Table_Final.rename(columns={Column:RepSeq}, inplace=True)
                continue

    # Sum the columns with the same name.
    Merged_Table_Final = Merged_Table_Final.groupby(axis=1, level=0).sum()

    if Transpose == True:
        Merged_Table_Final = Merged_Table_Final.T

    return Merged_Table_Final
    

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script merges OTU tables from DADA2 into a single one.\n'''
            '''It takes a representative fasta sequence and multiple OTU tables to be merged.\n'''
            '''The list of OTU tables should be separated by space, e.g., -l OTU_Table_1 OTU_Table_2 OTU_Table_3 ...\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [FastA File] -o [Output Table] -l [OTU_Table_1  OTU_Table_2...]\n'''
            '''Global mandatory parameters: -i [FastA File] -o [Output Table] -l [OTU Table List]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Input_Fasta', action='store', required=True, help='Input Fasta file with representative sequences')
    parser.add_argument('-l', '--list', dest='OTU_List', action='store', nargs='+', required=True, help='OTU table list, Separated by spaces.')
    parser.add_argument('-o', '--output', dest='Output_Table', action='store', required=True, help='Output OTU Table')
    parser.add_argument('--transpose', dest='Transpose', action='store_true', required=False, help='Output as sample in columns, OTU in rows, by default the other way.')
    args = parser.parse_args()

    Input_Fasta = args.Input_Fasta
    OTU_List = args.OTU_List
    Output_Table = args.Output_Table
    Transpose = args.Transpose

    Merged_OTU_Table = OTUTable_Merger(Input_Fasta, OTU_List, Transpose)

    Merged_OTU_Table.to_csv(Output_Table, sep="\t")

if __name__ == "__main__":
    main()