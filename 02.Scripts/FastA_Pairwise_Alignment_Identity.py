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
def child_initialize(_dictionary, _subsmatrix, _local):
     global sequence_dictionary, subsmatrix, local
     sequence_dictionary = _dictionary
     subsmatrix = _subsmatrix
     local = _local

def get_sequences(fasta_file):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    dictionary = {}
    with open(fasta_file) as Fasta:
        for title, sequence in SimpleFastaParser(Fasta):
            dictionary[title] = sequence
    return dictionary

def get_queries(query_file):
    query_list = []
    with open(query_file, 'r') as queries:
        for line in queries:
            query_list.append(line.strip())
    return query_list

def get_substitution_matrix():
    script_path = Path(__file__)
    script_dir = script_path.parent
    main_folder = Path(script_dir).parent
    martix_loc = main_folder / '00.Libraries/02.DNAfull_Sub_Matrix.txt'
    with open(martix_loc, 'rb') as filehandle:
        dnafull = pickle.load(filehandle)
    return dnafull


def perform_global_alignment(sequence_id):
    identity_dictionary = {}
    for title_b, sequence_b in sequence_dictionary.items():
        sequence_a = sequence_dictionary[sequence_id]
        alignment = pairwise2.align.globalds(sequence_a, sequence_b, subsmatrix,
        -10, -0.5, penalize_end_gaps=False, one_alignment_only=True)[0]
        if (sequence_id,title_b) in identity_dictionary or (title_b,sequence_id) in identity_dictionary:
            continue
        else:
            if local == False:
                identity = calculate_global_identity(alignment)
                identity_dictionary[(sequence_id,title_b)] = identity
            else:
                identity = calculate_local_identity(alignment)
                identity_dictionary[(sequence_id,title_b)] = identity
    outfile = sequence_id + '.id.txt'
    with open(outfile, 'w') as output:
        for pair, identity in identity_dictionary.items():
            output.write("{}\t{}\t{}\n".format(pair[0], pair[1], identity))
    print("Finished {}.".format(sequence_id))
    
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
    parser.add_argument('-f', '--fasta', dest='fasta_file', action='store', required=True, help='Fasta file with all sequences to align (references).')
    parser.add_argument('-q', '--query', dest='query', action='store', required=True, help='File with list of ids to use as queries vs all refs.')
    parser.add_argument('-t', '--threads', dest='threads', action='store', type=int, required=False, default=1, help='Threads to use. By default 1')
    parser.add_argument('--local', dest='local', action='store_true', required=False, help='Calculate local identities. By default calculates global.')
    args = parser.parse_args()

    fasta_file = args.fasta_file
    query = args.query
    local = args.local
    threads = args.threads

    subsmatrix = get_substitution_matrix()
    sequences = get_sequences(fasta_file)
    query_list = get_queries(query)
    try:
        pool = multiprocessing.Pool(threads, initializer = child_initialize, 
        initargs = (sequences, subsmatrix, local))
        pool.map(perform_global_alignment, query_list)
    finally:
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
