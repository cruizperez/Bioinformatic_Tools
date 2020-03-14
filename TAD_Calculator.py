#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.9
# Date:		   14 March 2020

# Description: This script parses a base-by-base sequencing depth file
# and calculates the TAD (Truncated Average Sequencing Depth) per genome.
# By default it calculates the TAD80 removing the 10% top and bottom
# covered bases, which takes care of highly (conserved) or poorly (contig
# edges) covered genome regions.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
from Bio.SeqIO.FastaIO import SimpleFastaParser
from statistics import mean
import argparse, sys

################################################################################

"""---2.0 Define Functions---"""
def calculate_tad_from_file(input_table, tad_percent, separator, outfile):
    to_remove = (100 - tad_percent)/2
    genome_seq = {}
    if separator is None:
        with open(input_table, 'r') as seqdepth:
            for line in seqdepth:
                if line.startswith('Sequence'):
                    continue
                else:
                    line = line.strip().split()
                    if line[0] not in genome_seq:
                        genome_seq[line[0]] = [line[2]]
                    else:
                        genome_seq[line[0]].append(line[2])
    else:
        with open(input_table, 'r') as seqdepth:
            for line in seqdepth:
                if line.startswith('Sequence'):
                    continue
                else:
                    line = line.strip().split()
                    genome = line[0].split(separator)[0]
                    if genome not in genome_seq:
                        genome_seq[genome] = [line[2]]
                    else:
                        genome_seq[genome].append(line[2])
    
    with open(outfile, 'w') as output:
        output.write("Genome\tTAD{}\n".format(tad_percent))
        for genome, depth in genome_seq.items():
            positions = round(len(depth)*to_remove/100)
            seq_sorted = sorted(depth, key=float)[positions:-positions]
            output.write("{}\t{}\n".format(genome, round(mean(seq_sorted),1)))

def calculate_tad_from_dict(input_dict, tad_percent, separator, outfile):
    to_remove = (100 - tad_percent)/2
    genome_seq = {}
    if separator is None:
        for contig, seq_depth in input_dict.items():
            if contig not in genome_seq:
                genome_seq[contig] = seq_depth
            else:
                genome_seq[contig] += seq_depth
    else:
        for contig, seq_depth in input_dict.items():
            genome = contig.split(separator)[0]
            if genome not in genome_seq:
                genome_seq[genome] = seq_depth
            else:
                genome_seq[genome] += seq_depth

    with open(outfile, 'w') as output:
        output.write("Genome\tTAD{}\n".format(tad_percent))
        for genome, depth in genome_seq.items():
            positions = round(len(depth)*to_remove/100)
            seq_sorted = sorted(depth, key=float)[positions:-positions]
            output.write("{}\t{}\n".format(genome, round(mean(seq_sorted),1)))


################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a base-by-base sequencing depth file\n'''
                        '''and calculates the TAD (Truncated Average Sequencing Depth) per genome.\n'''
                        '''By default it calculates the TAD80 removing the 10% top and bottom\n'''
                        '''covered bases, which takes care of highly (conserved) or poorly (contig\n'''
                        '''edges) covered genome regions.\n'''
                        '''Global mandatory parameters: [MagicBlast File] [Reference FastA]\n'''
                        'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input_seqdepth", dest='seqdepth', action='store', 
                        required=True, help="Input table with sequencing depth per base")
    parser.add_argument('-o', '--output_table', dest='output_table', action='store', 
                        required=True, help='Output table in the form [Sequence Name]\t[Position]\t[Depth]')
    parser.add_argument('--tad', dest='tad', action='store', required=False, default=80, type=int,
                        help='''TAD percentage to calculate. By default 80.\n''')
    parser.add_argument('--separator', dest='separator', action='store', required=False,
                        help='''String separating genome name from contig. By default None, i.e., calculates the sequencing depth per contig.\n'''
                             '''If separator has "-" or "--" pass it as --separator="--"''')
    args = parser.parse_args()

    seqdepth = args.seqdepth
    output_table = args.output_table
    tad = args.tad
    separator = args.separator

    # Calculate TAD and store results
    calculate_tad_from_file(seqdepth, tad, separator, output_table)


if __name__ == "__main__":
    main()
