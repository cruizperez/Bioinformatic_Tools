#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Institution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 14 August 2019

# Description: This script filters a MagicBlast tabular output based
# on the id percentage of the read and the length aligned.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Magic_Blast_Filter(Input_File, Alignment_Percent = 80, ID_Percent = 1):
    from pandas import DataFrame
    from numpy.random import randint

    Magic_Blast_Hits = {}
    with open(Input_File) as Tabular:
        for line in Tabular:
            line = line.strip().split()
            if '#' in line:
                continue
            elif float(line[2]) < ID_Percent:
                continue
            elif (float(line[7]) - float(line[6]) * 100 / float(line[15])) < Alignment_Percent:
                continue
            else:
                if line[0] not in Magic_Blast_Hits:
                    Magic_Blast_Hits[line[0]] = line[1:]
                else:
                    if float(line[12]) < float(Magic_Blast_Hits[line[0]][11]):
                        continue
                    elif float(line[12]) > float(Magic_Blast_Hits[line[0]][11]):
                        Magic_Blast_Hits[line[0]] = line[1:]
                    else:
                        if randint(2) > 0:
                            Magic_Blast_Hits[line[0]] = line[1:]
                        else:
                            continue
    Magic_Blast_DF = DataFrame.from_dict(Magic_Blast_Hits, orient='index')
    return Magic_Blast_DF
                    

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse
    from sys import argv
    import pandas as pd

    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script filters a MagicBlast tabular output based\n'''
            '''on the id percentage of the read and the length aligned.\n'''
            '''Usage: ''' + argv[0] + ''' -i [Input File] -o [Output File] -p [Identity] -f [Fraction Aligned]\n'''
            '''Global mandatory parameters: -i [Input File] -o [Output File] -p [Identity]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    parser.add_argument('-i', '--inputTab', dest='Input_Tab', action='store', required=True, help='Input MagicBLast tabular file.')
    parser.add_argument('-o', '--output', dest='Output', action='store', required=True, help='Output Table')
    parser.add_argument('-p', '--pidentity', dest='Percent_ID', action='store', required=True, type=int, help='Percentage identity of matches to retain.')
    parser.add_argument('-f', '--fraction', dest='Aln_Fraction', action='store', required=False, type=int, help='Minimum percentage of read aligned to be included. By default 80', default=80)
    args = parser.parse_args()

    Input_Tab = args.Input_Tab
    Output = args.Output
    Percent_ID = args.Percent_ID
    Aln_Fraction = args.Aln_Fraction

    MagicBLAST_Result = Magic_Blast_Filter(Input_Tab, Alignment_Percent = Aln_Fraction, ID_Percent = Percent_ID)
    MagicBLAST_Result.to_csv(Output, sep="\t", header=False)

if __name__ == "__main__":
    main()