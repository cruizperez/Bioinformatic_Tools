#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      0.9
# Date:         March 07, 2020

# Description: Calculates k-mer frequencies of every sequence in a given
# FastA file. Returns a matrix with sequences in rows and kmer in the columns.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
from itertools import product
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

################################################################################
"""---1.0 Define Functions---"""
def calculate_kmer_frequency(input_sequence_file, kmer_len):
    # Calculate possible nucleotide combinations
    nucl = "ATCG"
    possible_kmers = []
    possible_kmers_iter = product(nucl, repeat=kmer_len)
    for kmer in possible_kmers_iter:
        possible_kmers.append(''.join(kmer))
    # Retrieve sequences
    sequences = []
    with open(input_sequence_file) as fasta_input:
        for line in fasta_input:
            if line.startswith('>'):
                line = line.strip().split()[0].replace('>', '')
                sequences.append(line)
    # Create matrix
    kmer_frequency_matrix = pd.DataFrame(0, index=sequences, columns=possible_kmers)
    # Populate frequency matrix
    with open(input_sequence_file) as fasta_input:
        for identifier, sequence in SimpleFastaParser(fasta_input):
            kmers_sequence = get_kmers(sequence, kmer_len)
            for kmer in kmers_sequence:
                kmer_frequency_matrix.loc[identifier, kmer] += 1
    
    return kmer_frequency_matrix

def get_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
    return kmers


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script builds a sqlite database from a tab separated table.\n'''
            '''By default it assumes the first line of the input table has headers, if not\n'''
            '''you must provide a list of headers as --header header1,header2,header3...\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input Table] -d [Database Name]\n'''
            '''Global mandatory parameters: -i [Input Table] -d [Database Name]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_table', action='store', required=True,
                        help='Input tab-delimited table to parse, by default assumes headers are present')
    # parser.add_argument('-d', '--database', dest='database', action='store', required=True,
    #                     help='Database name, can be exisiting or new.')
    args = parser.parse_args()

    input_table = args.input_table
    # database = args.database

    # ----------------------------
    calculate_kmer_frequency(input_table, 4)
    # ----------------------------

if __name__ == "__main__":
    main()
