#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Institution:   Georgia Institute of Technology
# Version:	  0.9.0
# Date:		 29 July 2019

# Description: This script converts a FastA file into a FastQ file by
# assigning random qualities chosen by the user as lower and upper qualities.
# If qualities are not specified the highest quality is added to every base.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def FastA_to_FastQ(FastA_File, FastQ_File, LowQuality, HighQuality):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from random import choice
    # List of optional qualities.
    Qualities = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI'
    OutputFile = open(FastQ_File, 'w')
    with open(FastA_File) as FastA:
        for title, seq in SimpleFastaParser(FastA):
            QualScores = ''.join(choice(Qualities[LowQuality:HighQuality+1]) for i in range(len(seq)))
            OutputFile.write("@%s\n%s\n+\n%s\n" % (title, seq, QualScores))
    OutputFile.close()

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='This script converts a FastA file into a FastQ file by\n'
                                                    'assigning random qualities chosen by the user as lower and upper qualities.\n'
                                                    'If qualities are not specified the highest quality is added to every base.\n'
                                    '\nGlobal mandatory parameters: -i [Input FastA File]\n'
                                    '''\nOptional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='FastA_File', action='store', required=True, help='FastA file to convert')
    parser.add_argument('-o', '--output', dest='FastQ_File', action='store', required=False, help='Output FastQ, if not set "Sequences.fastq.', default="Sequences.fastq")
    parser.add_argument('--lower', dest='Low_Qual', action='store', required=False, help='Low quality value.', type=int, default=40)
    parser.add_argument('--higher', dest='High_Qual', action='store', required=False, help='High quality value.', type=int, default=40)

    args = parser.parse_args()

    FastA_File = args.FastA_File
    FastQ_File = args.FastQ_File
    Low_Qual = args.Low_Qual
    High_Qual = args.High_Qual

    # Run converter
    FastA_to_FastQ(FastA_File, FastQ_File, Low_Qual, High_Qual)

if __name__ == "__main__":
    main()
