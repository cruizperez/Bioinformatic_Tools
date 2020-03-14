#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.9
# Date:		   21 February 2020

# Description: This script calculates the global or local identities
# for each pair of sequences in a FastA-formatted alignment.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""

import sys, argparse, os
import multiprocessing

################################################################################

"""---2.0 Define Functions---"""
def child_initialize(_dictionary, _output):
     global sequence_dictionary, output
     sequence_dictionary = _dictionary
     output = _output

def Get_Alignment_Sequences(Alignment_File):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    Dictionary = {}
    with open(Alignment_File) as Fasta:
        for title, sequence in SimpleFastaParser(Fasta):
            Dictionary[title] = sequence
    return Dictionary

def calculate_global_identity(Sequence_ID):
    with open(Output, 'a') as Identity_File:
        for Title_B, Seq_B in Sequence_Dictionary.items():
            Seq_A = Sequence_Dictionary[Sequence_ID]
            Aln_Len = 0
            Match_Len = 0
            Len_Iter = max(len(Sequence_Dictionary[Sequence_ID]), len(Seq_B))
            for i in range(Len_Iter):
                if Seq_A[i] == "-" and Seq_B[i] == "-":
                    continue
                elif Seq_A[i] == Seq_B[i]:
                    Aln_Len += 1
                    Match_Len += 1
                else:
                    Aln_Len += 1
            Identity_File.write("{}\t{}\t{}\n".format(Sequence_ID,Title_B,Match_Len/Aln_Len))

def calculate_local_identity(Sequence_ID):
    with open(Output, 'a') as Identity_File:
        for Title_B, Seq_B in Sequence_Dictionary.items():
            Seq_A = Sequence_Dictionary[Sequence_ID]
            Aln_Len = 0
            Match_Len = 0
            Len_Iter = max(len(Sequence_Dictionary[Sequence_ID]), len(Seq_B))
            for i in range(Len_Iter):
                if Seq_A[i] == "-" and Seq_B[i] == "-":
                    continue
                elif Seq_A[i] == "-" and Seq_B[i] != "-":
                    continue
                elif Seq_A[i] != "-" and Seq_B[i] == "-":
                    continue
                elif Seq_A[i] == Seq_B[i]:
                    Aln_Len += 1
                    Match_Len += 1
                else:
                    Aln_Len += 1
            Identity_File.write("{}\t{}\t{}\n".format(Sequence_ID,Title_B,Match_Len/Aln_Len))

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script calculates the global or local identities\n'''
                        '''for each pair of sequences in a FastA-formatted alignment.\n'''
                        '''Global mandatory parameters: -f [Folder] -o [Output File] -i OR -l [Input files]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input_alignment', dest='input_aln', action='store', required=True, help='Alignment file in FastA format.')
    parser.add_argument('-o', '--output_file', dest='output_file', action='store', required=True, help='Tabular file to save identities.')
    parser.add_argument('-t', '--threads', dest='threads', action='store', type=int, required=False, default=1, help='Threads to use. By default 1')
    parser.add_argument('--local', dest='local', action='store_true', required=False, help='Calculate local identities. By default calculates global.')
    args = parser.parse_args()

    input_aln = args.input_aln
    output_file = args.output_file
    local = args.local
    threads = args.threads

    Sequences = Get_Alignment_Sequences(input_aln)
    Sequence_IDs = Sequences.keys()

    if local == True:
        try:
            pool = multiprocessing.Pool(threads, initializer = child_initialize, initargs = (Sequences, output_file))
            pool.map(calculate_local_identity, Sequence_IDs)
        finally:
            pool.close()
            pool.join()
    else:
        try:
            pool = multiprocessing.Pool(threads, initializer = child_initialize, initargs = (Sequences, output_file))
            pool.map(calculate_global_identity, Sequence_IDs)
        finally:
            pool.close()
            pool.join()

if __name__ == "__main__":
    main()
