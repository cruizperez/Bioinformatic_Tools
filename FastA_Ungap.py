#!/usr/bin/env python

################################################################################
"""---1.0 Import Modules---"""

import argparse, sys
from Bio.Seq import Seq
from Bio import SeqIO


################################################################################
"""---2.0 Define Functions---"""

def FastA_Ungapper(Fasta_File, Output_File):
    with open(Output_File, "w") as o:
        for record in SeqIO.parse(Fasta_File, "fasta"):
            record.seq = record.seq.ungap("-")
            SeqIO.write(record, o, "fasta")


################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''Removes gaps from an aligned FastA file'''
                                    'Global mandatory parameters: [FastA_File] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Output_File = args.Output_File

    FastA_Ungapper(Fasta_File, Output_File)

if __name__ == "__main__":
    main()
