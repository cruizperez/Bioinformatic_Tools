#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 24 March 2019

# Description: This script removed duplicate names and sequences from FastA files.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""
from random import randrange
import argparse, sys
from contextlib import ExitStack
from Bio.SeqIO.FastaIO import SimpleFastaParser

################################################################################
"""---2.0 Define Functions---"""

def FastA_Remove_Duplicate(Fasta_File, Output_File):
    Sequences = {}
    Output = open(Output_File, 'w')
    with open(Fasta_File) as Input:
        for title, seq in SimpleFastaParser(Input):
            if title in Sequences:
                pass
            else:
                Sequences[title] = seq
                Output.write(">%s\n%s\n" % (title, seq))
    Output.close()

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Remove duplicate sequences from a FastA file, by name'''
                                    'Global mandatory parameters: [Fasta File] [Output File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--inputFile", dest='Fasta_Files', required=True, help="Input FastA file")
    parser.add_argument("-o", "--outputFile", dest='Output_Files', required=True, help="Output FastA file")
    args = parser.parse_args()

    Fasta_Files = args.Fasta_Files
    Output_Files = args.Output_Files

    # Run remove duplicates Function
    FastA_Remove_Duplicate(Fasta_Files, Output_Files)

if __name__ == "__main__":
    main()
