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
def calculate_kmer_frequency(input_sequence_file, kmer_len, normalize):
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
        if normalize == True:
            for identifier, sequence in SimpleFastaParser(fasta_input):
                sequence_len = len(sequence)
                kmers_sequence = get_kmers(sequence, kmer_len)
                for kmer in kmers_sequence:
                    kmer_frequency_matrix.loc[identifier, kmer] += 1/sequence_len
        else:
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
            description='''Calculates k-mer frequencies of every sequence in a given FastA file.\n'''
                        '''Returns a matrix with sequences in rows and kmer in the columns.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [FastA File] -o [Output Table] -k [kmer lenght]\n'''
            '''Global mandatory parameters: -i [FastA File] -o [Output Table]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_fasta', action='store', required=True,
                        help='Input FastA file to parse.')
    parser.add_argument('-k', '--kmer', dest='kmer', action='store', required=False, type=int, default=4,
                        help='Kmer length, by default 4')
    parser.add_argument('-o', '--outfile', dest='output_table', action='store', required=True,
                        help='Output table')
    parser.add_argument('--normalize', dest='normalize', action='store_false', required=False,
                        help='Normalize frequencies by sequence length. By default True.')
    args = parser.parse_args()

    input_fasta = args.input_fasta
    kmer = args.kmer
    output_table = args.output_table
    normalize = args.normalize

    # ----------------------------
    kmer_frequency_table = calculate_kmer_frequency(input_fasta, kmer, normalize)
    kmer_frequency_table.to_csv(output_table, sep="\t", header=True, index=True)
    # ----------------------------

if __name__ == "__main__":
    main()
