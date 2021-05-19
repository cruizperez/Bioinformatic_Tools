#!/usr/bin/env python

"""
# ==============================================================================
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   April 30, 2021

# Description: 
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
def dereplicate_list(input_file, output_file):
    pairs_seen = {}
    with open(input_file, "r") as infile, open (output_file, "w") as outfile:
        for line in infile:
            line = line.strip().split()
            if (line[0], line[1]) in pairs_seen:
                continue
            elif (line[1], line[0]) in pairs_seen:
                continue
            else:
                outfile.write(f"{line[0]}\t{line[1]}\t{line[2]}\n")
                pairs_seen[(line[0], line[1])] = None
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, description=(
            f"This script remomes duplicated entries from a list with " 
            f"two columns and a third value."))
    mandatory_arguments = parser.add_argument_group("Mandatory arguments")
    mandatory_arguments.add_argument(
        '-i', '--input_file', dest='input_file', action='store', required=True,
        type=Path,
        help='Input file.')
    mandatory_arguments.add_argument(
        '-o', '--output', dest='output', action='store', required=True,
        help='Output file to store results.')
    # If no arguments were given
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Initialize input variables
    input_file = arguments.input_file
    output = arguments.output

    dereplicate_list(input_file, output)
# ==============================================================================


# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================