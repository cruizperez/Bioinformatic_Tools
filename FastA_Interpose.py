#!/usr/bin/env python

################################################################################
"""---1.0 Import Modules---"""

import argparse, sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
from contextlib import ExitStack

################################################################################
"""---2.0 Define Functions---"""

def FastA_Merger(FastaList, Output_File):
    Output_FH = open(Output_File, "w")
    with ExitStack() as stack:
        files = [stack.enter_context(open(fname)) for fname in FastaList]

        for seqs1, seqs2 in zip(SimpleFastaParser(files[0]), SimpleFastaParser(files[1])):
            Output_FH.write(">%s\n%s\n>%s\n%s" % (seqs1[0], seqs1[1], seqs2[0], seqs2[1]))

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Interposes two FastA files, usually paired end reads as Read_1.1, Read_1.2'''
                                    'Global mandatory parameters: [Cluster File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--inputFiles", dest='Input_Files', required=True, nargs=2, help="Input FastA paired end files")
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output paired end FastA file')
    args = parser.parse_args()

    Input_Files = args.Input_Files
    Output_File = args.Output_File

    # Run interpose Function
    FastA_Merger(Input_Files, Output_File)

if __name__ == "__main__":
    main()
