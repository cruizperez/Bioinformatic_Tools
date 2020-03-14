#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 24 March 2019

# Description: This script splits a FastA file into a fiven number of files.
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

def FastA_Splitter_Few(Fasta_File, Output_List):
    Num_Seqs = 0
    with open(Fasta_File) as Input:
        for line in Input:
            if line.startswith(">"):
                Num_Seqs += 1
    if Num_Seqs < len(Output_List):
        Output_List = Output_List[0:Num_Seqs]
    with ExitStack() as stack:
        files = [stack.enter_context(open(fname, 'w')) for fname in Output_List]
        Counter = len(files)
        File_num = 0
        with open(Fasta_File) as Input:
            for title, seq in SimpleFastaParser(Input):
                if File_num == Counter:
                    File_num = 1
                    files[File_num-1].write(">%s\n%s\n" % (title, seq))
                else:
                    File_num += 1
                    files[File_num-1].write(">%s\n%s\n" % (title, seq))

def FastA_Splitter_Several(Fasta_File, Output_List):
    Num_Seqs = 0
    with open(Fasta_File) as Input:
        for line in Input:
            if line.startswith(">"):
                Num_Seqs += 1
    if Num_Seqs < len(Output_List):
        Output_List = Output_List[0:Num_Seqs]
    if (Num_Seqs - int(Num_Seqs/len(Output_List)) * len(Output_List)) != 0:
        Seqs_per_File = int(Num_Seqs/len(Output_List)) + 1
    else:
        Seqs_per_File = int(Num_Seqs/len(Output_List))

    with open(Fasta_File) as Input:
        Index = 0
        Count = 0
        for title, seq in SimpleFastaParser(Input):
            if Count < (Seqs_per_File):
                CurrentFile = open(Output_List[Index], 'w+')
                CurrentFile.write(">%s\n%s\n" % (title, seq))
                CurrentFile.close()
                Count += 1
            else:
                Index += 1
                CurrentFile = open(Output_List[Index], 'w+')
                CurrentFile.write(">%s\n%s\n" % (title, seq))
                CurrentFile.close()
                Count = 1


################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Splits FastA file into a given number of files'''
                                    'Global mandatory parameters: [Fasta File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--inputFiles", dest='Fasta_Files', required=True, help="Input FastA file")
    parser.add_argument('-p', '--prefix', dest='Prefix', action='store', help='Output prefix, by default "Fasta_File_"')
    parser.add_argument('-n', '--number', dest='Number', action='store', type=int, default=2, help='Number of files to split into, 2 by default')
    args = parser.parse_args()

    Fasta_Files = args.Fasta_Files
    Prefix = args.Prefix
    Number = args.Number

    if Prefix == None:
        Prefix = "Fasta_File_"
    # Run split Function
    Output_list = [Prefix + str(i) + ".fa" for i in range(1, Number+1)]
    if Number <= 1000:
        FastA_Splitter_Few(Fasta_Files, Output_list)
    else:
        FastA_Splitter_Several(Fasta_Files, Output_list)

if __name__ == "__main__":
    main()
