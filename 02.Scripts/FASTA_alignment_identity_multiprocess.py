#!/usr/bin/env python

"""
# ==============================================================================
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   April 30, 2021

# Description: This script calculates the global or local identities
# of all pairwise alignments in a FASTA file. The identity is defined as the
# number of matches over the total alignment length (global) or over the
# nucleotides present (excluding gaps; local).
# ==============================================================================
"""


# ==============================================================================
# Import modules
# ==============================================================================
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path
from sys import argv
from sys import exit

import multiprocessing
import argparse
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
def child_initialize(_sequence_dictionary, _target_sequences):
     global sequence_dictionary, target_sequences
     sequence_dictionary = _sequence_dictionary
     target_sequences = _target_sequences

def get_alignment_sequences(Alignment_File):
    Dictionary = {}
    with open(Alignment_File) as Fasta:
        for title, sequence in SimpleFastaParser(Fasta):
            title = title.strip().split()[0]
            Dictionary[title] = sequence
    return Dictionary

# def calculate_global_identity(Sequence_ID):
#     with open(output, 'a') as Identity_File:
#         for Title_B, Seq_B in sequence_dictionary.items():
#             Seq_A = sequence_dictionary[Sequence_ID]
#             Aln_Len = 0
#             Match_Len = 0
#             Len_Iter = max(len(sequence_dictionary[Sequence_ID]), len(Seq_B))
#             for i in range(Len_Iter):
#                 if Seq_A[i] == "-" and Seq_B[i] == "-":
#                     continue
#                 elif Seq_A[i] == Seq_B[i]:
#                     Aln_Len += 1
#                     Match_Len += 1
#                 else:
#                     Aln_Len += 1
#             Identity_File.write("{}\t{}\t{}\n".format(Sequence_ID,Title_B,Match_Len/Aln_Len))

def calculate_local_identity(query_id):
    identity_list = []
    query_seq = sequence_dictionary[query_id]
    for target_id in target_sequences:
        target_seq = sequence_dictionary[target_id]
        alignment_lenght = 0
        match_length = 0
        length_longest = max(len(query_seq), len(target_seq))
        for i in range(length_longest):
            if query_seq[i] == "-" or target_seq[i] == "-":
                continue
            elif query_seq[i] == target_seq[i]:
                alignment_lenght += 1
                match_length += 1
            else:
                alignment_lenght += 1
        if alignment_lenght > 0:
            identity = match_length/alignment_lenght
        else:
            identity = 'na'
        identity_list.append((query_id, target_id, identity))
    return identity_list
# ==============================================================================



# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, description=(
            f"This script calculates the global or local identity for each " 
            f"pairwise alignment in a FASTA file.\nIdentities are defined as " 
            f"the number of matches over the total alignment length or over " 
            f"the nucleotides present (excluding gaps).\n" 
            f"Mandatory parameters: -f [input FASTA] -o [output file]\n" 
            f"Optional parameters: See '{argv[0]} -h'\n"))
    mandatory_arguments = parser.add_argument_group("Mandatory arguments")
    mandatory_arguments.add_argument(
        '-f', '--fasta', dest='fasta_file', action='store', required=True,
        type=Path,
        help='FASTA file with all sequences to align.')
    mandatory_arguments.add_argument(
        '-o', '--output', dest='output', action='store', required=True,
        help='Output file to store identity results.')
    optional_arguments = parser.add_argument_group("Optional arguments")
    optional_arguments.add_argument(
        '--ql', dest='query_list', action='store', required=False,
        nargs='*',
        help=(
            f'Space-separated list of sequences to use as queries.'
            f'By default it uses all sequences.'))
    optional_arguments.add_argument(
        '--qf', dest='query_file', action='store', required=False,
        help=(
            f'File with list of sequences to use as queries (one per line). '
            f'By default it uses all sequences. Has priority over "--ql".'))
    optional_arguments.add_argument(
        '--tl', dest='target_list', action='store', required=False,
        nargs='*',
        help=(
            f'Space-separated list of sequences to use as targets.'
            f'By default it uses all sequences.'))
    optional_arguments.add_argument(
        '--tf', dest='target_file', action='store', required=False,
        help=(
            f'File with list of sequences to use as targets (one per line). '
            f'By default it uses all sequences. Has priority over "--tl".'))
    optional_arguments.add_argument(
        '-t', '--threads', dest='threads',
        action='store', type=int, required=False, default=1,
        help='Threads to use. By default 1')
    optional_arguments.add_argument(
        '--local', dest='local', action='store_true', required=False,
        help='Calculate local identities. By default calculates global.')
    # If no arguments were given
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()
    # Initialize input variables
    fasta_file = arguments.fasta_file
    output = arguments.output
    query_list = arguments.query_list
    query_file = arguments.query_file
    target_list = arguments.target_list
    target_file = arguments.target_file
    threads = arguments.threads
    local = arguments.local

    fasta_sequences = get_alignment_sequences(fasta_file)
    if query_list == None:
        query_list = fasta_sequences.keys()
    if target_list == None:
        target_list = fasta_sequences.keys()

    if local:
        try:
            pool = multiprocessing.Pool(
                threads, initializer = child_initialize, initargs = (
                    fasta_sequences, target_list))
            alignment_identities = pool.map(
                calculate_local_identity, query_list)
        finally:
            pool.close()
            pool.join()
    # else:
    #     try:
    #         pool = multiprocessing.Pool(threads, initializer = child_initialize, initargs = (Sequences, output_file))
    #         alignment_identities = pool.map(calculate_global_identity, Sequence_IDs)
    #     finally:
    #         pool.close()
    #         pool.join()

    # print(alignment_identities)
    with open(output, 'w') as outfile:
        for result in alignment_identities:
            for identity in result:
                outfile.write(f"{identity[0]}\t{identity[1]}\t{identity[2]}\n")

if __name__ == "__main__":
    main()
