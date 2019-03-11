#!/usr/bin/env python

################################################################################
"""---1.0 Import Modules---"""

import argparse, sys
from Bio.SeqIO.FastaIO import SimpleFastaParser


################################################################################
"""---2.0 Define Functions---"""

def Coordinates_Parser(Coordinates_File):
    Coordinates = {}
    with open(Coordinates_File) as Input:
        for line in Input:
            line = line.strip().split(sep="\t")
            # Find and sort the positions of the sequences.
            Positions = sorted(list(map(int, line[1:3])))
            # Add positions to dictionary and also add strand location.
            Coordinates[line[0]] = Positions
            Coordinates[line[0]].append(line[3])
    return Coordinates


def FastA_Sequence_Extract(Fasta_File, Coordinates_Dictionary, Output_File):
    Output_FH = open(Output_File, "w")
    with open(Fasta_File) as Input:
        for title, seq in SimpleFastaParser(Input):
            if title in Coordinates_Dictionary:
                # Extract sequence from contig in the given positions from the dictionary.
                sequence = Seq(seq[int(Coordinates_Dictionary[title][0])-1:int(Coordinates_Dictionary[title][1])-1])
                # Evaluate location on strand.
                if Coordinates_Dictionary[title][2] == "-":
                    Output_FH.write(">%s\n%s\n" % (title, secuencia.reverse_complement()))
                    print(title, sequence.reverse_complement())
                else:
                    print(title, secuencia)
                    Output_FH.write(">%s\n%s\n" % (title, sequence))

################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''Extracts sequence(s) from a FastA file based on the positions along the sequence.
                                                    The positions file should be in the form [Sequence ID] [Pos1] [Pos2], and can also have information
                                                    on the strand (as + or -) in which case, the reverse complement is returned for those in the complementary strand'''
                                    'Global mandatory parameters: [FastA_File] [Coordinates_File] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    parser.add_argument('-c', '--coord', dest='Coordinates_File', action='store', required=True, help='Coordinates file')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Output_File = args.Output_File
    Coordinates_File = args.Coordinates_File

    Coordinates = Coordinates_Parser(Coordinates_File)
    FastA_Sequence_Extract(Fasta_File, Coordinates, Output_File)

if __name__ == "__main__":
    main()
