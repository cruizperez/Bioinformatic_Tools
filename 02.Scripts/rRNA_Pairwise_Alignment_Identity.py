#!/usr/bin/env python

"""
# ==============================================================================
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   April 30, 2021

# Description: This script performs a global pairwise alignment using the 
# Needleman–Wunsch algorithm implemented in Biopython and using the same
# parameters as in EMBOSS_Needle alignment (DNAfull matrix).
# It then calculates the global or local identities, defined as the number
# of matches over the total alignment length or over the nucleotides present
# (excluding gaps).
# ==============================================================================
"""


# ==============================================================================
# Import modules
# ==============================================================================
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tempfile import TemporaryDirectory
from Bio import pairwise2
from pathlib import Path
from typing import Tuple
from typing import Dict
from typing import List
from typing import Any
from sys import argv
from sys import exit

import numpy as np
import argparse
import random
import pickle
import string
import shutil
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
def child_initialize(
    _dictionary: Dict[str,str], _subsmatrix: Dict[str,np.ndarray],
    _target_sequences: List[str], _local: bool,
    _query_as_target: bool, _temp_directory: str) -> None:
     global sequence_dictionary, subsmatrix, local, target_sequences
     global query_as_target, temp_directory
     sequence_dictionary = _dictionary
     subsmatrix = _subsmatrix
     target_sequences = _target_sequences
     local = _local
     query_as_target = _query_as_target
     temp_directory = _temp_directory


def get_sequences(fasta_file: Path) -> Dict[str,str]:
    dictionary = {}
    with open(fasta_file) as Fasta:
        for title, sequence in SimpleFastaParser(Fasta):
            title = title.strip().split()[0]
            dictionary[title] = sequence
    return dictionary


def get_substitution_matrix() -> Dict[str,np.ndarray]:
    script_path = Path(__file__)
    script_dir = script_path.parent
    main_folder = Path(script_dir).parent
    martix_loc = main_folder / '00.Libraries/02.DNAfull_Sub_Matrix.txt'
    with open(martix_loc, 'rb') as filehandle:
        dnafull = pickle.load(filehandle)
    return dnafull


def perform_global_alignment(
    queries: List[str], targets: List[str], rrna_sequences: Dict[str,str],
    subsmatrix: Any, local: bool, output: Path) -> None:
    with open(output, "w") as outfile:
        for query in queries:
            print(f"Processing {query}")
            query_seq = rrna_sequences[query]
            for target in targets:
                target_seq = rrna_sequences[target]
                # Perform pairwise alignments
                alignment = pairwise2.align.globalds(
                    query_seq, target_seq, subsmatrix,
                    -10, -0.5, penalize_end_gaps=False,
                    one_alignment_only=True)[0]
                if local:
                    identity = calculate_local_identity(alignment)
                else:
                    identity = calculate_global_identity(alignment)
                outfile.write(f"{query}\t{target}\t{identity}\n")
                outfile.flush()
    # query_seq = sequence_dictionary[query]
    # random_name = random_string_generator(10)
    # output = f"{temp_directory.name}/{query}_{random_name}.temp"
    # print(query, "\n", query_seq)
    # with open(output, "w") as outfile:
    #     for target in target_sequences:
    #         target_seq = sequence_dictionary[target]
    #         if 'nnn' in target_seq.lower() or \
    #             'nnn' in query_seq.lower():
    #             continue
    #         # Perform pairwise alignments
    #         alignment = pairwise2.align.globalds(
    #             query_seq, target_seq, subsmatrix,
    #             -10, -0.5, penalize_end_gaps=False,
    #             one_alignment_only=True)[0]
    #         # Calculate identities
    #         if local:
    #             identity = calculate_global_identity(alignment)
    #         else:
    #             identity = calculate_local_identity(alignment)
    #         if query_as_target:
    #             outfile.write(f"{target}\t{query}\t{identity}\n")
    #         else:
    #             outfile.write(f"{query}\t{target}\t{identity}\n")
    # return output


def calculate_global_identity(alignment):
    aln_len = alignment[4]
    ident_bases = 0
    seq_a = alignment[0]
    seq_b = alignment[1]
    for i in range(aln_len):
        if seq_a[i] == "-" and seq_b[i] == "-":
            continue
        elif seq_a[i] == seq_b[i]:
            ident_bases += 1
        else:
            continue
    return round(ident_bases/aln_len, 3)


def calculate_local_identity(alignment):
    glob_aln_len = alignment[4]
    ident_bases = 0
    matched_aln_len = 0
    seq_a = alignment[0]
    seq_b = alignment[1]
    for i in range(glob_aln_len):
        if seq_a[i] == "-" or seq_b[i] == "-":
            continue
        elif seq_a[i] == seq_b[i]:
            ident_bases += 1
            matched_aln_len += 1
        else:
            matched_aln_len += 1
    return round(ident_bases/matched_aln_len, 3)


def random_string_generator(str_size: int):
    chars = string.ascii_letters + string.digits
    return ''.join(random.choice(chars) for x in range(str_size))
 

# ==============================================================================


# ==============================================================================
# Define main functions
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, description=(
            f"This script performs a global pairwise alignment using the"
            f"Needleman–Wunsch algorithm implemented in Biopython.\n"
            f"It uses the same parameters as in EMBOSS_Needle alignment "
            f"(DNAfull matrix).\nIt then calculates the global or local" 
            f"identity for each pairwise alignment.\nIdentities are defined as " 
            f"the number of matches over the total alignment length or over " 
            f"the nucleotides present (excluding gaps).\n" 
            f"Mandatory parameters: -f [input FASTA] -o [output file]\n" 
            f"Optional parameters: See '{argv[0]} -h'\n"))
    mandatory_arguments = parser.add_argument_group("Mandatory arguments")
    mandatory_arguments.add_argument(
        '-f', '--fasta', dest='fasta_file', action='store', required=True,
        type=Path,
        help='FASTA file with all sequences to align.')
    mandatory_arguments.add_argument(
        '-o', '--output', dest='output', action='store', required=True,
        help='Output file to store identity results.')
    optional_arguments = parser.add_argument_group("Optional arguments")
    optional_arguments.add_argument(
        '--ql', dest='query_list', action='store', required=False,
        nargs='*',
        help=(
            f'Space-separated list of sequences to use as queries.'
            f'By default it uses all sequences.'))
    optional_arguments.add_argument(
        '--qf', dest='query_file', action='store', required=False,
        help=(
            f'File with list of sequences to use as queries (one per line). '
            f'By default it uses all sequences. Has priority over "--ql".'))
    optional_arguments.add_argument(
        '--tl', dest='target_list', action='store', required=False,
        nargs='*',
        help=(
            f'Space-separated list of sequences to use as targets.'
            f'By default it uses all sequences.'))
    optional_arguments.add_argument(
        '--tf', dest='target_file', action='store', required=False,
        help=(
            f'File with list of sequences to use as targets (one per line). '
            f'By default it uses all sequences. Has priority over "--tl".'))
    optional_arguments.add_argument(
        '-t', '--threads', dest='threads',
        action='store', type=int, required=False, default=1,
        help='Threads to use. By default 1')
    optional_arguments.add_argument(
        '--local', dest='local', action='store_true', required=False,
        help='Calculate local identities. By default calculates global.')
    # If no arguments were given
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()
    # Initialize input variables
    fasta_file = arguments.fasta_file
    output = arguments.output
    query_list = arguments.query_list
    query_file = arguments.query_file
    target_list = arguments.target_list
    target_file = arguments.target_file
    threads = arguments.threads
    local = arguments.local

    # ==========================================================================
    # Process and validate input variables
    # ==========================================================================
    # Validate input fasta
    if not fasta_file.is_file():
        exit(
            f"ERROR: Input FASTA file cannot be found. "
            f"Check the location and try again.")
    # Validate query lists
    query_sequences = None
    if (query_list and query_file) or query_file:
        query_sequences = []
        with open(query_file, "r") as infile:
            for line in infile:
                query_sequences.append(line.strip())
    elif query_list:
        query_sequences = query_list
    # Validate target lists
    target_sequences = None
    if (target_list and target_file) or target_file:
        target_sequences = []
        with open(target_file, "r") as infile:
            for line in infile:
                target_sequences.append(line.strip())
    elif target_list:
        target_sequences = target_list
    # ==========================================================================


    # ==========================================================================
    # Process and validate input variables
    # ==========================================================================
    # Import substitution matrix
    subsmatrix = get_substitution_matrix()
    # Parse sequences from FASTA file
    sequences = get_sequences(fasta_file)
    # If not queries or targets are provided get them from the whole set
    if not query_sequences:
        query_sequences = list(sequences.keys())
    if not target_sequences:
        target_sequences = list(sequences.keys())
    # # Get longer sequence to take advantage of multiprocessing
    # if len(target_sequences) < len(query_sequences):
    #     sequences_global = target_sequences
    #     sequences_multiprocess = query_sequences
    #     query_as_target = False
    # elif len(target_sequences) > len(query_sequences):
    #     sequences_global = query_sequences
    #     sequences_multiprocess = target_sequences
    #     query_as_target = True
    # else:
    #     sequences_global = target_sequences
    #     sequences_multiprocess = query_sequences
    #     query_as_target = False
    # ==========================================================================
    perform_global_alignment(
    queries=query_sequences, targets=target_sequences, rrna_sequences=sequences,
    subsmatrix=subsmatrix, local=local, output=output)

    # ==========================================================================
    # Create temporary directory to store results
    # ==========================================================================
    # temp_directory = TemporaryDirectory(dir="./")
    # ==========================================================================


    # ==========================================================================
    # Initialize global shared variables and perform alignments
    # ==========================================================================
    # try:
    #     multiprocessing.set_start_method("fork")
    #     pool = multiprocessing.Pool(threads, initializer = child_initialize, 
    #     initargs = (
    #         sequences, subsmatrix, sequences_global, local, query_as_target,
    #         temp_directory))
    #     alignment_results = pool.map(
    #         perform_global_alignment, sequences_multiprocess)
    #     alignment_results.wait()
    # finally:
    #     pool.close()
    #     pool.join()
    
    # print("Merging results")
    # with open(output, 'w') as outfile:
    #     for file in alignment_results:
    #         with open(file) as temporary:
    #             shutil.copyfileobj(temporary, outfile)
# ==============================================================================


# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================