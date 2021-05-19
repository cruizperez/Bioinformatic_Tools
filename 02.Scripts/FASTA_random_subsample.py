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
from functools import partial
from typing import List
from typing import Tuple
from random import sample
from pathlib import Path
from sys import argv
from sys import exit

import multiprocessing
import argparse
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Function to read FASTA and get list of sequence ids
def get_sequence_ids(fasta_file: Path) -> Tuple[List[str], Path]:
    """
    Gets identifiers for all sequences in a FASTA file

    Args:
        fasta_file (Path): Input FASTA file

    Returns:
        List[str]: List of sequence identifiers
    """
    sequence_list = []
    with open(fasta_file) as infile:
        for line in infile:
            if line.startswith(">"):
                line = line.strip().split()[0]
                identifier = line.replace(">", "")
                sequence_list.append(identifier)
    return sequence_list, fasta_file

# Function to get a random subsample from FASTA file
def get_random_subsample(
    sequence_lists: Tuple[List[str]], percentage: int) -> None:
    """
    Randomnly subsamples FASTA file to a given percentage of the total
    initial sequences

    Args:
        fasta_file (Path): Input FASTA file
        output (Path): Output subsampled FASTA file
        input_list (List[str]): List of initial sequence identifiers
        percentage (int): Target percentage to subsample
    """
    # Get total number of sequences
    input_list = sequence_lists[0] # Get initial list of ids
    fasta_file = Path(sequence_lists[1]) # Get path of original fasta file
    output =  Path(fasta_file.name).with_suffix(f".{percentage}.fasta")
    total_seqs = len(input_list)
    num_seqs = int(total_seqs * percentage / 100)
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
            f"Mandatory parameters: -i [input FASTA] -p [percetage]"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory arguments")
    mandatory_arguments.add_argument(
        "-i", "--input_fasta", dest='input_fasta', action='store', nargs='+',
        required=True, help="Input FASTA file")
    mandatory_arguments.add_argument(
        '-p', '--percentage', dest='percentage', action='store', type=int, 
        required=True, help='Target percentage to subsample')
    # Setup optional arguments
    optional_arguments = parser.add_argument_group("Optional arguments")
    optional_arguments.add_argument(
        "-t", "--threads", dest='threads', action='store', type=int, 
        required=False, default=1, help="Threads to use. Default 1.")
    # If no arguments were given
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    input_fasta = arguments.input_fasta
    percentage = arguments.percentage
    threads = arguments.threads

    # Get total sequence list
    try:
        pool = multiprocessing.Pool(threads)
        sequence_lists = pool.map(get_sequence_ids, input_fasta)
    finally:
        pool.close()
        pool.join()
    # Subsample fasta files
    try:
        pool = multiprocessing.Pool(threads)
        pool.map(partial(
            get_random_subsample, percentage=percentage), sequence_lists)
    finally:
        pool.close()
        pool.join()
# ==============================================================================


# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
