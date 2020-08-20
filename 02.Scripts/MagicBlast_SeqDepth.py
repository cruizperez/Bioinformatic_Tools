#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   1.0
# Date:		   26 January 2020

# Description: This script parses a MagicBlast tabular output and
# the reference sequences (in FastA format) and returns the base by base
# sequencing depth of each reference contig/genome. 
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse, sys

################################################################################

"""---2.0 Define Functions---"""
def get_genome_sizes(fasta_file):
    """
    Calculates the length of each sequence within a FastA file
    
    Arguments:
        fasta_file {filepath} -- Fasta file with reference sequences
    
    Returns:
        [dictionary] -- Lengths per sequence
    """
    genome_sizes = {}
    with open(fasta_file) as fasta_input:
        for title, seq in SimpleFastaParser(fasta_input):
            genome_sizes[title] = len(seq)
    return genome_sizes


def calculate_seq_depth(magicblast_file, genome_sizes):
    """
    Calculates the base-by-base sequencing depth.
    
    Arguments:
        magicblast_file {filepath} -- File with MagicBlast tabular output.
        genome_sizes {dictionary} -- Lengths per sequence
    
    Returns:
        [dictionary] -- Sequencing depth per DNA fragment
    """
    genome_seq_depth = {}
    
    with open(magicblast_file) as magicblast:
        for line in magicblast:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split()
                sequence = line[1]
                seq_start = min(int(line[8]), int(line[9]))
                seq_end = max(int(line[8]), int(line[9]))
                if line[1] not in genome_seq_depth:
                    genome_seq_depth[sequence] = np.zeros(genome_sizes[sequence])
                    genome_seq_depth[sequence][seq_start-1:seq_end-1] += 1
                else:
                    genome_seq_depth[sequence][seq_start-1:seq_end-1] += 1
    
    return  genome_seq_depth

def calculate_seq_depth_sorted(magicblast_file, genome_sizes, output_table):
    """
    Calculates the base-by-base sequencing depth from a sorted file for lower memory consumption.
    
    Arguments:
        magicblast_file {filepath} -- File with MagicBlast tabular output.
        genome_sizes {dictionary} -- Lengths per sequence
        output_table {filepath} -- Output table file
    """
    current_contig = None
    current_bases = None

    with open(output_table, 'w') as output, open(magicblast_file, 'r') as magicblast:
        output.write("Sequence\tPosition\tDepth\n")
        for line in magicblast:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split()
                sequence = line[1]
                seq_start = min(int(line[8]), int(line[9]))
                seq_end = max(int(line[8]), int(line[9]))
                if current_contig == None:
                    current_contig = sequence
                    current_bases = np.zeros(genome_sizes[sequence])
                    current_bases[seq_start-1:seq_end-1] += 1
                elif current_contig == sequence:
                    current_bases[seq_start-1:seq_end-1] += 1
                elif current_contig != sequence:
                    for base, coverage in np.ndenumerate(current_bases):
                        output.write("{}\t{}\t{}\n".format(current_contig, str(base[0]+1), str(coverage)))
                    current_contig = sequence
                    current_bases = np.zeros(genome_sizes[sequence])
                    current_bases[seq_start-1:seq_end-1] += 1

#TODO: Fix this table saving to avoid creating two huge data structures (save line by line).
def save_sequencing_depth_table(genome_seq_depth, output_table):
    """
    Saves a dictionary of arrays as a table with
    Sequence_Name   Position    Depth
    
    Arguments:
        genome_seq_depth {dictionary} -- Array of per base sequence depth per genome
        output_table {path} -- Path of file to save table
    """
    sequences = []
    positions = []
    depths = []
    genome_seq_depth_table = pd.DataFrame(columns=["Sequence","Position","Depth"])
    for sequence, depth_array in genome_seq_depth.items():
        sequence_list = [sequence] * len(depth_array)
        sequences += sequence_list
        position_array = np.array(range(len(depth_array))) + 1
        positions += list(position_array)
        depths += list(depth_array)
    genome_seq_depth_table["Sequence"] = sequences
    genome_seq_depth_table["Position"] = positions
    genome_seq_depth_table["Depth"] = depths
    genome_seq_depth_table.to_csv(output_table, sep="\t", header=True, index=False)


################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a MagicBlast tabular output and\n'''
                        '''the reference sequences (in FastA format) and returns\n'''
                        '''the base by base sequencing depth of each reference contig/genome\n'''
                        'Global mandatory parameters: -i [MagicBlast File] -f [Reference FastA] -o [Output Table]\n'
                        'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input_magicblast", dest='magic_blast', action='store', 
                        required=True, help="Input MagiBlast tabular output")
    parser.add_argument('-f', '--fasta_sequences', dest='fasta_sequences', action='store', 
                        required=True, help='FastA file of reference sequences')
    parser.add_argument('-o', '--output_table', dest='output_table', action='store', 
                        required=True, help='Output table in the form [Sequence Name]\t[Position]\t[Depth]')
    parser.add_argument('-s', '--sorted', dest='sorted_input', action='store_true', 
                        required=False, help='If input is sorted by the second column this will save memory.')
    args = parser.parse_args()

    magic_blast = args.magic_blast
    fasta_sequences = args.fasta_sequences
    output_table = args.output_table
    sorted_input = args.sorted_input

    # Calculate Genome Length and Sequencing Depth
    genome_sizes = get_genome_sizes(fasta_sequences)
    if sorted_input == True:
        calculate_seq_depth_sorted(magic_blast, genome_sizes, output_table)
    else:
        genome_seq_depth = calculate_seq_depth(magic_blast, genome_sizes)
        # Save output
        save_sequencing_depth_table(genome_seq_depth, output_table)

if __name__ == "__main__":
    main()
