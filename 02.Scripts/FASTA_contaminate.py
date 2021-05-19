#!/usr/bin/env python

"""
# ==============================================================================
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   April 30, 2021

# Description: This script takes a FASTA file a source of proteins and 
# randomly adds a given number of chosen proteins to another FASTA file
# ==============================================================================
"""


# ==============================================================================
# Import modules
from Bio.SeqIO.FastaIO import SimpleFastaParser
from typing import Dict
from random import sample
from pathlib import Path
from sys import argv
from sys import exit

import argparse
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Function to read source FASTA and get dictionary of sequences
def get_source_sequences(fasta_file: Path) -> Dict[str, str]:
    """
    Get source sequences from a FASTA file and stores in a dictionary

    Args:
        fasta_file (Path): Input FASTA file

    Returns:
        Dict[str, str]: Dictionary with source sequences
    """
    source_sequences = {}
    with open(fasta_file) as infile:
        for title, sequence in SimpleFastaParser(infile):
            source_sequences[title] = sequence
    return source_sequences

# Function to get a random subsample from FASTA file
def add_proteins_to_fasta(
    target_fasta: Path, source_sequences: Dict[str, str],
    proteins_to_add: int, output: Path) -> None:
    # Get list of possible proteins to add
    source_ids = source_sequences.keys()
    # Randomly select proteins
    selected_ids = sample(source_ids, proteins_to_add)
    # Write final FASTA file
    with open(target_fasta, 'r') as infile, open(output, 'w') as outfile:
        for line in infile:
            outfile.write(f"{line}")
        for seq_id in selected_ids:
            outfile.write(f">{seq_id}_dupl\n{source_sequences[seq_id]}\n")
# ==============================================================================         


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script randomly adds sequences to a FASTA file "
            f"from a source FASTA file.\n"
            f"Mandatory parameters: -i [input FASTA] -s [source FASTA] "
            f"-o [output FASTA] -p [proteins to add]"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory arguments")
    mandatory_arguments.add_argument(
        "-i", "--input_fasta", dest='input_fasta', action='store',
        required=True, help="Input FASTA file")
    mandatory_arguments.add_argument(
        '-s', '--source', dest='source', action='store', 
        required=True, help='FASTA with proteins to be added (source)')
    mandatory_arguments.add_argument(
        '-o', '--output', dest='output', action='store', 
        required=True, help='Output FASTA file')
    mandatory_arguments.add_argument(
        '-p', '--proteins', dest='proteins', action='store', type=int, 
        required=True, help='Number of proteins to add')
    # If no arguments were given
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    input_fasta = Path(arguments.input_fasta)
    source = Path(arguments.source)
    output = Path(arguments.output)
    proteins = arguments.proteins

    # Get initial source proteins
    source_proteins = get_source_sequences(source)
    # Add proteins to original FASTA file
    add_proteins_to_fasta(
        input_fasta, source_proteins, proteins, output)
# ==============================================================================


# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
