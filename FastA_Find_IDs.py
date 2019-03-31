#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import sys, argparse, os
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

"""----------------------------- 1.0 Define Functions -----------------------------"""

def FastA_Filter_List(List, FastaFile, Output):
    Seq_ID_list = []
    with open(FastaFile) as Fasta_Input:
        for title, seq in SimpleFastaParser(Fasta_Input):
            Seq_ID_list.append(title.split()[0])
    Output_List = open(Output, 'w')
    with open(List) as Seq_IDs:
        for line in Seq_IDs:
            line = line.strip()
            if line[0] in Seq_ID_list:
                Output_List.write("%s\tYes\n" % (line[0]))
            else:
                Output_List.write("%s\tNo\n" % (line[0]))
    Output_File.close()


### ------------------------------- Main function ------------------------------

def main():
    parser = argparse.ArgumentParser(description='''Parses a file with sequence IDs and tells if they are present in a FastA file'''
                                    'Global mandatory parameters: [FastA_File] [Output_File] [ID List File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    parser.add_argument('-l', '--list', dest='ID_List', action='store', required=True, help='List of IDs to filter')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Output_File = args.Output_File
    ID_List = args.ID_List

    FastA_Filter_List(ID_List, Fasta_File, Output_File)

if __name__ == "__main__":
    main()
