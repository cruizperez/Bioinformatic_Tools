#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import sys, argparse, os
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

"""----------------------------- 1.0 Define Functions -----------------------------"""

def FastA_Filter_List(List, FastaFile, Reverse, Output):
    Seq_ID_list = []
    Records = 0
    with open(List) as Seq_IDs:
        for line in Seq_IDs:
            line = line.strip()
            Seq_ID_list.append(line)
            Records += 1
    Fasta_out = open(Output, 'w')
    with open(FastaFile) as Fasta_in:
        if Reverse == True:
            print("Excluding " + str(Records) + " records from output")
            for title, seq in SimpleFastaParser(Fasta_in):
                if title.split()[0] not in Seq_ID_list:
                    Fasta_out.write(">%s\n%s\n" % (title, seq))
        else:
            print("Retrieving " + str(Records) + " records from input")
            for title, seq in SimpleFastaParser(Fasta_in):
                if title.split()[0] in Seq_ID_list:
                    Fasta_out.write(">%s\n%s\n" % (title, seq))
    Fasta_out.close()


def main():
    parser = argparse.ArgumentParser(description='''Filter a FastA file based ona provided list of IDs. It can exclude or retrieve the sequences using the --reverse flag'''
                                    'Global mandatory parameters: [FastA_File] [Output_File] [ID List File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    parser.add_argument('-l', '--list', dest='ID_List', action='store', required=True, help='List of IDs to filter')
    parser.add_argument('--reverse', action='store_true', help='Exclude the sequences in the list file. By default False, i.e. retrieves those in the list')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Output_File = args.Output_File
    ID_List = args.ID_List
    Reverse = args.reverse

    FastA_Filter_List(ID_List, Fasta_File, Reverse, Output_File)

if __name__ == "__main__":
    main()
