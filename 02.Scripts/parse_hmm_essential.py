#!/usr/bin/env python

"""
# ==============================================================================
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   April 30, 2021

# Description: This script randomly subsamples a FASTA file to a given
# percentage of the original sequences. 
# ==============================================================================
"""


# ==============================================================================
# Import modules
from typing import List
from random import sample
from pathlib import Path
from sys import argv
from sys import exit
import argparse
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Function to read FASTA and get list of sequence ids
def get_sequence_ids(fasta_file: Path) -> List[str]:
    sequence_list = []
    with open(fasta_file) as infile:
        for line in infile:
            if line.startswith(">"):
                line = line.strip().split()[0]
                identifier = line.replace(">", "")
                sequence_list.append(identifier)
    return sequence_list

# Function to get a random subsample from FASTA file
def get_random_subsample(
    fasta_file: Path, output: Path,
    input_list: List[str], percentage: int) -> None:
    # Get total number of sequences
    total_seqs = len(input_list)
    num_seqs = total_seqs * percentage / 100
    selected_ids = sample(input_list, num_seqs)
    with open(fasta_file, "r") as infile, open(output, "w") as outfile:
        retain = False
        for line in infile:
            if line.startswith(">"):
                retain = False
                identifier = line.strip().split()[0]
                identifier = identifier.replace(">", "")
                if identifier in selected_ids:
                    retain = True
                    outfile.write(f"{line.strip()}\n")
            else:
                if retain:
                    outfile.write(f"{line.strip()}\n")
                else: 
                    continue
# ==============================================================================         


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script randomly subsamples a FASTA file "
            f"to a given percentage.\n"
            f"Mandatory parameters: -i [input FASTA] -o [output FASTA] "
            f"-p [percetage]"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory arguments")
    mandatory_arguments.add_argument(
        "-i", "--input_fasta", dest='input_fasta', action='store', 
        required=True, help="Input FASTA file")
    mandatory_arguments.add_argument(
        '-o', '--output_fasta', dest='output_fasta', action='store', 
        required=True, help='Output FASTA to store results')
    mandatory_arguments.add_argument(
        '-p', '--percentage', dest='percentage', action='store', type='int', 
        required=True, help='Target percentage to subsample')
    # If no arguments were given
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    input_fasta = arguments.input_fasta
    output_fasta = arguments.output_fasta
    percentage = arguments.percentage

    # Get total sequence list
    sequence_list = get_sequence_ids(input_fasta)
    get_random_subsample(input_fasta, output_fasta, sequence_list, percentage)
# ==============================================================================


# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
