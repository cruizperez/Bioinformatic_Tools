#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 Oct 04 2019

# Description: This script parses a CheckM "bin_stats_ext.tsv" file from
the analyze workflow along with a list of SCGs and returns a matrix with 
each genome and SCG found in it.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Parse_CheckM_Result(CheckM_Output, OutputMatrix):
    import pandas as pd
    import ast

    Genome_List = []
    Genome_Information = {}
    with open(CheckM_Output, 'r') as Input:
        for line in Input:
            line = line.strip().split("\t")
            Genome_List.append(line[0])
            Genome_Information[line[0]] = ast.literal_eval(line[1])
    Genome_List = list(set(Genome_List))

    Result_Matrix = pd.DataFrame(index=Genome_List)
    for Genomes, Information in Genome_Information.items():
        for copies in enumerate(['GCN0', 'GCN1', 'GCN2', 'GCN3', 'GCN4', 'GCN5+']):
            for SC in Information[copies[1]]:
                Result_Matrix.loc[Genomes, SC] = copies[0]

    Result_Matrix.to_csv(OutputMatrix, sep="\t")


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse
    from sys import argv

    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a CheckM "bin_stats_ext.tsv" file from\n'''
                        '''the analyze workflow along with a list of SCGs and returns a matrix with\n'''
                        '''each genome and SCG found in it.\n'''
            '''Usage: ''' + argv[0] + ''' -i [CheckM Output] -l [SCG list] -o [Gene Copy Matrix]\n'''
            '''Global mandatory parameters: -i [CheckM Output] -l [SCG list] -o [Gene Copy Matrix]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='CheckM_Output', action='store', required=True, help='Input CheckM Result.')
    parser.add_argument('-o', '--output', dest='OutputMatrix', action='store', required=True, help='Output Matrix')
    args = parser.parse_args()

    CheckM_Output = args.CheckM_Output
    OutputMatrix = args.OutputMatrix

    Parse_CheckM_Result(CheckM_Output, OutputMatrix)

if __name__ == "__main__":
    main()