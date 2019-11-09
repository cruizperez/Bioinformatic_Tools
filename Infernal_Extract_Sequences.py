#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 October 31 2019

# Description: This script extracts the RNA sequences of a genome/contig
# predicted using Infernal from its tabular output.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Sequence_Extract(Input, Genome_Folder, Output):
    from pathlib import Path
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.Seq import Seq
    import sys

    Filename = Path(Genome_Folder)
    rRNA_Location = {}
    # Read infernal file and extract RNA locations.
    with open(Input, 'r') as Infernal:
        File = None
        Type = None
        Contig = None
        Start = None
        End = None
        Strand = None
        for line in Infernal:
            if '#' not in line:
                line = line.strip().split()
                Type = line[1]
                Contig = line[3]
                Start = int(line[9])
                End = int(line[10])
                Strand = line[11]
            elif 'Query file' in line:
                line = line.strip().split()
                File = line[3]
            # Add entry to dictionary as File = [Contig, Start, End, Strand]
        if Contig == None:
            sys.exit('No hits present in {}'.format(Input))
        if File not in rRNA_Location:
            rRNA_Location[File] =[[Contig, Start, End, Strand, Type]]
        else:
            rRNA_Location[File].append([Contig, Start, End, Strand, Type])
    # For each entry (file) open the corresponding fasta file and extract the sequence.
    Output_Handle = open(Output, 'w')
    for File, Matches in rRNA_Location.items():
        Counter =1
        Genome_Name = File.split(".")
        Genome_Name = '.'.join(Genome_Name[:-1])
        for Entry in Matches:
            New_Path = Filename / File
            with New_Path.open() as Fasta:
                for title, sequence in SimpleFastaParser(Fasta):
                    title = title.split()[0]
                    if title == Entry[0]:
                        Start = min(Entry[1], Entry[2])
                        End = max(Entry[1], Entry[2])
                        Sequence = Seq(sequence[Start-1:End-1])
                        if Entry[3] == "-":
                            Sequence = Sequence.reverse_complement()
            Output_Handle.write(">{}-{}_{}\n{}\n".format(Genome_Name, Entry[4], Counter, Sequence))
            Counter += 1
    Output_Handle.close()

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script extracts the RNA sequences of a genome/contig\n'''
            '''predicted using Infernal from its tabular output and outputs as a Fasta file\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Infernal Tab File] -o [Output FastA] -g [Genome Folder]\n'''
            '''Global mandatory parameters: -i [Infernal Tab File] -o [Output FastA] -g [Genome Folder]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--inf', dest='Infernal', action='store', required=True, help='Infernal tabular file')
    parser.add_argument('-o', '--out', dest='Output_Fasta', action='store', required=True, help='Output FastA file')
    parser.add_argument('-g', '--gen', dest='Genomes', action='store', required=True, help='Folder with genome fasta files')
    args = parser.parse_args()

    Infernal = args.Infernal
    Output_Fasta = args.Output_Fasta
    Genomes = args.Genomes

    Sequence_Extract(Infernal, Genomes, Output_Fasta)

if __name__ == "__main__":
    main()


