#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.9
# Date:		   21 February 2020

# Description: This script performs a global pairwise alignment using the 
# Needleman–Wunsch algorithm implemented in Biopython and using the same
# parameters as in EMBOSS_Needle alignment (DNAfull matrix).
# It then calculates the global or local identities, defined as the number
# of matches over the total alignment length or over the nucleotides present
# (excluding gaps).
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""

import sys, argparse, os
import multiprocessing
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pickle
from pathlib import Path

################################################################################

"""---2.0 Define Functions---"""
def child_initialize(_dictionary, _output, _subsmatrix, _local):
     global sequence_dictionary, output, subsmatrix, local
     sequence_dictionary = _dictionary
     output = _output
     subsmatrix = _subsmatrix
     local = _local

def get_sequences(fasta_file):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    dictionary = {}
    with open(fasta_file) as Fasta:
        for title, sequence in SimpleFastaParser(Fasta):
            dictionary[title] = sequence
    return dictionary

def get_substitution_matrix():
    script_path = Path(__file__)
    script_dir = script_path.parent
    main_folder = Path(script_dir).parent
    martix_loc = main_folder / '00.Libraries/02.DNAfull_Sub_Matrix.txt'
    with open(martix_loc, 'rb') as filehandle:
        dnafull = pickle.load(filehandle)
    return dnafull


def perform_global_alignment(sequence_id):
    for title_b, sequence_b in sequence_dictionary.items():
        sequence_a = sequence_dictionary[sequence_id]
        alignment = pairwise2.align.globalds(sequence_a, sequence_b, subsmatrix,
        -10, -0.5, penalize_end_gaps=False, one_alignment_only=True)[0]
        if local == False:
            identity = calculate_global_identity(alignment)
        else:
            identity = calculate_local_identity(alignment)
        with open(output, 'a') as outfile:
            outfile.write("{}\t{}\t{}\n".format(sequence_id,title_b,identity))
    
    
def calculate_global_identity(alignment):
    aln_len = alignment[4]
    ident_bases = 0
    seq_a = alignment[0]
    seq_b = alignment[1]
    for i in range(aln_len):
        if seq_a[i] == "-" and seq_b[i] == "-":
            continue
        elif seq_a[i] == seq_b[i]:
            ident_bases += 1
        else:
            continue
    return round(ident_bases/aln_len, 3)

def calculate_local_identity(alignment):
    glob_aln_len = alignment[4]
    ident_bases = 0
    matched_aln_len = 0
    seq_a = alignment[0]
    seq_b = alignment[1]
    for i in range(glob_aln_len):
        if seq_a[i] == "-" or seq_b[i] == "-":
            continue
        elif seq_a[i] == seq_b[i]:
            ident_bases += 1
            matched_aln_len += 1
        else:
            matched_aln_len += 1
    return round(ident_bases/matched_aln_len, 3)

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script performs a global pairwise alignment using the\n''' 
                        '''Needleman–Wunsch algorithm implemented in Biopython and using the same\n'''
                        '''parameters as in EMBOSS_Needle alignment (DNAfull matrix).\n'''
                        '''It then calculates the global or local identities, defined as the number\n'''
                        '''of matches over the total alignment length or over the nucleotides present\n'''
                        '''(excluding gaps).\n'''
                        '''Global mandatory parameters: -i [Input Fasta] -o [Output File]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_file', action='store', required=True, help='Fasta file with sequences to align.')
    parser.add_argument('-o', '--output_file', dest='output_file', action='store', required=True, help='Tabular file to save identities.')
    parser.add_argument('-t', '--threads', dest='threads', action='store', type=int, required=False, default=1, help='Threads to use. By default 1')
    parser.add_argument('--local', dest='local', action='store_true', required=False, help='Calculate local identities. By default calculates global.')
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    local = args.local
    threads = args.threads

    subsmatrix = get_substitution_matrix()
    sequences = get_sequences(input_file)
    sequence_ids = sequences.keys()
    try:
        pool = multiprocessing.Pool(threads, initializer = child_initialize, 
        initargs = (sequences, output_file, subsmatrix, local))
        pool.map(perform_global_alignment, sequence_ids)
    finally:
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
