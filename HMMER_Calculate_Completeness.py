#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   07 February 2020

# Description: This script calculates the completeness and contamination
# of genomes from a filtered HMM search results file. 
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
from random import randint

################################################################################

"""---2.0 Define Functions---"""
def calculate_completeness(hmmsearch_file, hmm_models, outfile, genome_separator):
    number_scg = 0
    with open(hmm_models, 'r') as models:
        for line in models:
            if "NAME" in line:
                number_scg += 1

    genome_completeness = {}
    with open(hmmsearch_file, 'r') as input_hmm:
        for line in input_hmm:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split()
                genome = line[0].split(genome_separator)[0]
                scg_accession = line[3]
                if genome not in genome_completeness:
                    genome_completeness[genome] = {scg_accession: 1}
                else:
                    if scg_accession not in genome_completeness[genome]:
                        genome_completeness[genome][scg_accession] = 1
                    else:
                        genome_completeness[genome][scg_accession] += 1
    
    with open(outfile, 'w') as output:
        output.write("Genome\tCompleteness\tContamination\n")
        for genomes, scgs in genome_completeness.items():
            redundant = 0
            for copies in scgs.values():
                if copies > 1:
                    redundant += 1
            output.write("{}\t{}\t{}\n".format(genomes, round((len(scgs)/number_scg)*100, 1), round(redundant/number_scg*100, 1)))

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a hmmsearch (HMMER) tabular ouput and\n'''
                        '''returns the best hit per gene based on the best domain OR sequence bitscore.\n'''
                        '''Global mandatory parameters: [HMMSearch Tab File] [Output File]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input_hmmsearch", dest='hmmsearch', action='store', 
                        required=True, help="Input hmmsearch tabular output")
    parser.add_argument('-m', '--models_file', dest='models_file', action='store', 
                        required=True, help='File to store filtered hmmsearch results')
    parser.add_argument('-o', '--outfile', dest='outfile', action='store', 
                        required=True, help='File to store filtered hmmsearch results')
    parser.add_argument('--genome_separator', dest='genome_separator', action='store', default="--", 
                        required=False, help='Filter hits by domain score, by default True')
    args = parser.parse_args()

    hmmsearch = args.hmmsearch
    models_file = args.models_file
    outfile = args.outfile
    genome_separator = args.genome_separator

    # Calculate completeness and contamination
    calculate_completeness(hmmsearch, models_file, outfile, genome_separator)

if __name__ == "__main__":
    main()
