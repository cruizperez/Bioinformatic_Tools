#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

"""----------------------------- 1.0 Define Functions -----------------------------"""

def FastA_Filter(FastaFile, Length, Output, max=False):
    with open(FastaFile) as Input:
        with open(Output) as Output:
            if max = True
                for title, seq in SimpleFastaParser(Input):
                    if len(seq) <= int(Length):
                        Ouput.write(">%s\n%s\n" % (title, seq))
            else:
                for title, seq in SimpleFastaParser(Input):
                    if len(seq) >= int(Length):
                        Ouput.write(">%s\n%s\n" % (title, seq))


def main():
    parser = argparse.ArgumentParser(description='''Filter a FastA file based on the minimum or maximum length required per sequence (by default min)'''
                                    'Global mandatory parameters: [FastA_File] [Output_File] [Length]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    parser.add_argument('-l', '--length', dest='Length', action='store', required=True, help='Minumum (or Maximum) length')
    parser.add_argument('--max', action='store_true', help='Retrieves sequences of less than max length, by default false')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Output_File = args.Output_File
    Length = args.Length
    Max = args.max

    FastA_Filter(Fasta_File, Length, Output_File, Max)

if __name__ == "__main__":
    main()
