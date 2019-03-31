#!/usr/bin/env python

################################################################################
"""---1.0 Import Modules---"""

import argparse, sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
from random import randrange


################################################################################
"""---2.0 Define Functions---"""

def FastA_Extract_Longest(Fasta_File, Output_File):
    Output_FH = open(Output_File, "w")
    Length = 0
    with open(Fasta_File) as Input:
        for title, seq in SimpleFastaParser(Input):
            if len(seq) < Length:
                continue
            elif len(seq) == Length:
                if randrange(0, 2) > 0:
                    ID_Name = title
                    sequence = seq
                else:
                    continue
            else:
                Length = len(seq)
                ID_Name = title
                sequence = seq
    Output_FH.write(">%s\n%s\n" % (ID_Name, sequence))

################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''Extracts the longest sequence in a FastA file. If two sequences are equally long, it selects one at random'''
                                    'Global mandatory parameters: [FastA_File] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Output_File = args.Output_File

    FastA_Extract_Longest(Fasta_File, Output_File)

if __name__ == "__main__":
    main()
