#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import sys, argparse, os
from Bio import SeqIO
import pandas as pd

"""----------------------------- 1.0 Define Functions -----------------------------"""

def HitConfidence(GeneID, identity, bitscore, long, qlen = 100, slen = 100):
    if long == True:
        if (int(min(qlen,slen))/int(max(qlen,slen))*100 >= 80) and float(bitscore) >= 80 and float(identity) >= 40:
            return(GeneID)
    else:
        if float(bitscore) >= 80 and float(identity) >= 40:
            return(GeneID)

def Blast_Parser(BlastFile):
    print("Reading Blast Output")
    BlastFile = pd.read_csv(BlastFile, sep="\t", header=None)
    if len(BlastFile.columns) > 12:
        print("I am assuming your Blast output has qlen and slen besides the standard output columns")
        ID_List = (BlastFile.apply(lambda row: HitConfidence(row[0], row[2], row[11], True, row[12], row[13]), axis=1)).tolist()
    else:
        print("I am assuming your Blast output has the standard output")
        ID_List = (BlastFile.apply(lambda row: HitConfidence(row[0], row[2], row[11], True), axis=1)).tolist()
    ID_List = [i for i in ID_List if i is not None]
    return(ID_List)

def FastA_Filter(List, FastaFile, Reverse, Output):
    records = SeqIO.parse(FastaFile, "fasta")
    if Reverse == True:
        with open(Output, 'a') as f_out:
            for record in records:
                if record.id not in List:
                    SeqIO.write(record, f_out, "fasta")
    if Reverse == False:
        with open(Output, 'a') as f_out:
            for record in records:
                if record.id in List:
                    SeqIO.write(record, f_out, "fasta")

def main():
    parser = argparse.ArgumentParser(description='''Given a Blast output and a FastA file, determines which sequences had good matches and retrieves
                                                    either those with good matches or those without for further processing'''
                                    'Global mandatory parameters: [FastA_File] [Blast_File] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-b', '--blast', dest='Blast_File', action='store', required=True, help='Blast output of the FastA file search agains a DB')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    parser.add_argument('--inverse', action='store_false', help='Retrieve the sequences with good matches. By default False, i.e. retrieves those without matches')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Blast_File = args.Blast_File
    Output_File = args.Output_File
    Inverse = args.inverse

    ID_List = Blast_Parser(Blast_File)
    FastA_Filter(ID_List, Fasta_File, Inverse, Output_File)

if __name__ == "__main__":
    main()
