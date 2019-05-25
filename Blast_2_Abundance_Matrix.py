#!/usr/bin/env python

################################################################################
"""---1.0 Import Modules---"""

import pandas as pd
import pathlib
import argparse, sys

################################################################################
"""---2.0 Define Functions---"""

def Blast_2_Matrix(Files_List, Num_Ext):
    Num_Ext = int(Num_Ext)
    # Create empty dataframe.
    Read_Matrix = pd.DataFrame()

    for File in Files_List:
        print("Processing {}...").format(File)
        Filename = pathlib.Path(File)
        # Remove as many extensions as indicated.
        for num in range(1,Num_Ext+1):
            Name = Filename.stem
            Filename = pathlib.Path(Name)
        # Add column and fill with 0.
        Read_Matrix[Filename] = 0
        with open(File) as Blast_File:
            for line in Blast_File:
                line = line.strip()
                line = line.split(sep="\t")
                if line[1] not in Read_Matrix.index:
                    Read_Matrix = Read_Matrix.reindex(Read_Matrix.index.values.tolist()+[line[1]], fill_value=0)
                    Read_Matrix.loc[line[1],Filename] += 1
                else:
                    Read_Matrix.loc[line[1],Filename] += 1

    return Read_Matrix

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Creates a read mapping table based on multiple blast outputs'''
                                    'Global mandatory parameters: [Blast Files (list)]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-l', '--List', dest='Files_List', action='store', required=True, nargs='+', help='List of Blast Files (Name will be used as column names)')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=False, help='Output table with reads mapped per subject, if none, "Matrix_Counts.tab".', default="Matrix_Counts.tab")
    parser.add_argument('--ext', dest='Num_Ext', action='store', type=int, help='Number of extensions to remove from file name, e.g. 1 to remove .blast from Genome1.blast. If none, 1', default=1)
    args = parser.parse_args()

    Files_List = args.Files_List
    Output_File = args.Output_File
    Num_Ext = args.Num_Ext

    # Return dataframes
    Read_Table = Blast_2_Matrix(Files_List, Num_Ext)

    # Export both dataframes to tab-separated tables.
    Read_Table.to_csv(Output_File, sep='\t')

if __name__ == "__main__":
    main()
