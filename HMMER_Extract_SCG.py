#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   10 February 2020

# Description: This script parses a hmmsearch (HMMER) tabular ouput and
# extracts the fasta sequences per single copy marker gene (SCG).
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
from Bio.SeqIO.FastaIO import SimpleFastaParser
import multiprocessing
from functools import partial

################################################################################

"""---2.0 Define Functions---"""
def hmm_extract_scg_genes(hmmsearch_file):

    scg_groups = {}
    with open(hmmsearch_file, 'r') as input:
        for line in input:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split()
                if line[3] in scg_groups:
                    scg_groups[line[3]].append(line[0])
                else:
                    scg_groups[line[3]] = [line[0]]
    return scg_groups

def scg_extract_sequence(scg_group, information):
    outfile = scg_group[0] + ".faa"
    proteins = scg_group[1]
    sequence_file = information[0]
    contig_separator = information[1]
    with open(sequence_file) as fasta_input, open(outfile, 'w') as fasta_output:
        for title, sequence in SimpleFastaParser(fasta_input):
            if title in proteins:
                genome_id = title.split(contig_separator)[0:-1]
                fasta_output.write(">{}\n{}\n".format(genome_id, sequence))
            else:
                continue

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a hmmsearch (HMMER) tabular ouput and\n'''
                        '''extracts the fasta sequences per single copy marker gene (SCG).\n'''
                        '''Global mandatory parameters: -i [HMMSearch Tab File] -f [Protein FastA]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input_hmmsearch", dest='hmmsearch', action='store', 
                        required=True, help="Input hmmsearch tabular output")
    parser.add_argument('-f', '--fasta_file', dest='fasta_file', action='store', 
                        required=True, help='File to store filtered hmmsearch results')
    parser.add_argument('-t', '--threads', dest='threads', action='store', type=int, 
                        required=False, help='Threads to use. By default 1')
    parser.add_argument('--separator', dest='separator', action='store', 
                        required=False, help='Contig delimiter. By default "--", e.g. Genome1--contig1_gene1')
    args = parser.parse_args()

    hmmsearch = args.hmmsearch
    fasta_file = args.fasta_file
    threads = args.threads
    separator = args.separator

    scg_groups = hmm_extract_scg_genes(hmmsearch)
    scg_list = []
    for key, value in scg_groups.items():
        scg_list.append((key, value))

    try:
        pool = multiprocessing.Pool(threads)
        pool.map(partial(scg_extract_sequence, information=(fasta_file, separator)), scg_list)
    finally:
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
