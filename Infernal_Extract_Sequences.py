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

def Sequence_Extract(Input, Output, Genome_File=None):
    from pathlib import Path
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.Seq import Seq
    import sys

    if Genome_File != None:
        Genome_File = Path(Genome_File)
    rRNA_Location = {}
    # Read infernal file and extract RNA locations.
    with open(Input, 'r') as Infernal:
        for line in Infernal:
            if '#' not in line:
                line = line.strip().split()
                Type = line[2]
                Contig = line[0]
                Start = int(line[7])
                End = int(line[8])
                Strand = line[9]
                # Add entry to dictionary as File = [Contig, Start, End, Strand]
                if Contig not in rRNA_Location:
                    rRNA_Location[Contig] =[[Contig, Start, End, Strand, Type]]
                else:
                    rRNA_Location[Contig].append([Contig, Start, End, Strand, Type])
            elif Genome_File == None and 'Target file' in line:
                line = line.strip().split()
                Genome_File = Path(line[3])
            
        if len(rRNA_Location) == 0:
            sys.exit('No hits present in {}'.format(Input))
    
    # For each entry (file) open the corresponding fasta file and extract the sequence.
    Output_Handle = open(Output, 'w')
    for Contig, Matches in rRNA_Location.items():
        Counter = 1
        for Entry in Matches:
            with Genome_File.open() as Fasta:
                for title, sequence in SimpleFastaParser(Fasta):
                    title = title.split()[0]
                    if title == Contig:
                        Start = min(Entry[1], Entry[2])
                        End = max(Entry[1], Entry[2])
                        Sequence = Seq(sequence[Start-1:End])
                        if Entry[3] == "-":
                            Sequence = Sequence.reverse_complement()
            Output_Handle.write(">{}-{}_{}\n{}\n".format(Contig, Entry[4], Counter, Sequence))
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
            '''Usage: ''' + sys.argv[0] + ''' -i [Infernal Tab File] -o [Output FastA] -g [Genome File]\n'''
            '''Global mandatory parameters: -i [Infernal Tab File] -o [Output FastA]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--inf', dest='Infernal', action='store', required=True, help='Infernal tabular file')
    parser.add_argument('-o', '--out', dest='Output_Fasta', action='store', required=True, help='Output FastA file')
    parser.add_argument('-g', '--gen', dest='Genome', action='store', required=False, help='Genome fasta file, by default inferred from input')
    args = parser.parse_args()

    Infernal = args.Infernal
    Output_Fasta = args.Output_Fasta
    Genome = args.Genome

    Sequence_Extract(Infernal, Output_Fasta, Genome)

if __name__ == "__main__":
    main()


