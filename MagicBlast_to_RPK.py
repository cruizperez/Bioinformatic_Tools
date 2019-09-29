#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Institution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 14 August 2019

# Description: This script filters a MagicBlast tabular output based
# on the id percentage of the read and the length aligned.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

################################################################################

"""---2.0 Define Functions---"""

def Get_Genome_Sizes(Fasta_File):
    Genome_Sizes = {}
    with open(Fasta_File) as Fasta:
        for title, seq in SimpleFastaParser(Fasta):
            if "VIRSorter" in title:
                Genome = title
            else:
                Genome = title.split("_")
                Genome = "_".join(Genome[:-1])
            if Genome not in Genome_Sizes:
                Genome_Sizes[Genome] = len(seq)
    return Genome_Sizes


def Calculate_Seq_Depth(MagicBlast_File, Genome_Sizes):
    Genomes = {}
    Read_Lenght = []
    
    with open(MagicBlast_File) as Input:
        for line in Input:
            line = line.strip().split()
            if "VIRSorter" in line[1]:
                Genome = line[1]
            else:
                Genome = line[1].split("_")
                Genome = "_".join(Genome[:-1])
            if Genome not in Genomes:
                Genomes[Genome] = 0
                Genomes[Genome] += 1
                Read_Lenght.append(int(line[15]))
            else:
                Genomes[Genome] += 1
                Read_Lenght.append(int(line[15]))
    
    Read_Lenght = np.asarray(Read_Lenght)
    Mean_Read_Len = np.mean(Read_Lenght)
    
    return(Genomes, Mean_Read_Len)

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''This script parses MAPLE's output and returns a table with each KEGG KO ID with its associated module and step in that module\n
                                                     The input should be any .matrix file from either the KAAS or the BLAST folders.'''
                                    'Global mandatory parameters: [Maple_Output] [Output Table]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--inputMagicBlast", dest='MagicBlast', action='store', required=True, help="Input Maple Annotation")
    parser.add_argument('-f', '--genomes', dest='Genomes_File', action='store', required=True, help='Mapping file from KO IDs to gene, downloaded from https://www.genome.jp/kegg-bin/get_htext#B1')
    parser.add_argument('-o', '--outputTable', dest='Output_Table', action='store', required=True, help='Output table in the form KEGG_ID Module Step')
    parser.add_argument('-s', '--sample', dest='Sample_Name', action='store', required=True, help='Output table in the form KEGG_ID Module Step')
    args = parser.parse_args()

    MagicBlast = args.MagicBlast
    Genomes_File = args.Genomes_File
    Output_Table = args.Output_Table
    Sample_Name = args.Sample_Name


    # Calculate Genome Length and Sequencing Depth
    Genome_Sizes = Get_Genome_Sizes(Genomes_File)
    (Seq_Depth, Mean_Read_Len) = Calculate_Seq_Depth(MagicBlast, Genome_Sizes)
    
    # Merge into a single Table
    Genome_Depth = {}
    for key, value in Seq_Depth.items():
        Genome_Depth[key] = ((value * Mean_Read_Len) / (Genome_Sizes[key])) * 1000


    Table = pd.DataFrame.from_dict(Genome_Depth, orient='index')
    Table.columns = [Sample_Name]
    Table.index.name = "Genome"
    Table.to_csv(Output_Table, sep="\t")


if __name__ == "__main__":
    main()
