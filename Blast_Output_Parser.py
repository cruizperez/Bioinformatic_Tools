#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import sys, argparse, os
from Bio import SeqIO
import pandas as pd

"""----------------------------- 1.0 Define Functions -----------------------------"""
### ------------------------------Match filter--------------------------------------
def HitConfidence(line, id, bitscore, evalue, Aln_Percent = None):
    if Aln_Percent == None:
        if float(line[2]) >= float(id) and float(line[11]) >= float(bitscore) and float(line[10]) <= float(evalue):
            High_Quality_Match = True
        else:
            High_Quality_Match = False
    else:
        if (int(line[3])*100/int(line[12]) >= Aln_Percent) and float(line[11]) >= float(bitscore) and float(line[2]) >= float(id) and float(line[10]) <= float(evalue):
            High_Quality_Match = True
        else:
            High_Quality_Match = False
    return(High_Quality_Match)

### ----------------------------- Blast Parser ----------------------------------------
def Blast_Parser(BlastFile, Output, id, bitscore, evalue, Aln_Percent = None):
    Blast_Dict = {}
    print("Reading Blast Output")
    # Check if Blast output has qlen and slen in addition to std output.
    with open(BlastFile) as BlastFile_Input:
        if len(BlastFile_Input.readline().strip().split("\t")) == 14:
            print("I am assuming your Blast output has qlen and slen besides the standard output columns")
            long = True
        else:
            print("I am assuming your Blast output has the standard output")
            long = False
    # Check if any parameter was given for match filtering.
    if all(variable is None for variable in [id, bitscore, evalue, Aln_Percent]):
        print("Only retrieving best match per sequence.")
        # Give only best match based on bitscore.
        with open(BlastFile) as BlastFile_Input:
                for line in BlastFile_Input:
                    line = line.strip().split("\t")
                    if line[0] not in Blast_Dict:
                        Blast_Dict[line[0]] = line[1:-1]
                    else:
                        if float(line[11]) >= float(Blast_Dict[line[0]][10]):
                            Blast_Dict[line[0]] = line[1:-1]
                        elif float(line[11]) == float(Blast_Dict[line[0]][10]):
                                if randrange(0, 2) > 0:
                                    Blast_Dict[line[0]] = line[1:-1]
                        else:
                            pass
    else:
        print("Performing match filtering and retrieving best match based on the parameters you provided.")
        # Do match filtering based on parameters provided and retrieve best match based on best bitscore.
        with open(BlastFile) as BlastFile_Input:
            if long == True:
                for line in BlastFile_Input:
                    line = line.strip().split("\t")
                    Good = HitConfidence(line, id, evalue, bitscore, Aln_Percent)
                    if Good == True:
                        if line[0] not in Blast_Dict:
                            Blast_Dict[line[0]] = line[1:-1]
                        else:
                            if float(line[11]) >= float(Blast_Dict[line[0]][10]):
                                Blast_Dict[line[0]] = line[1:-1]
                            elif float(line[11]) == float(Blast_Dict[line[0]][10]):
                                if randrange(0, 2) > 0:
                                    Blast_Dict[line[0]] = line[1:-1]
                            else:
                                pass
                    else:
                        pass
            else:
                for line in BlastFile_Input:
                    line = line.strip().split("\t")
                    Good = HitConfidence(line, id, evalue, bitscore)
                    if Good == True:
                        if line[0] not in Blast_Dict:
                            Blast_Dict[line[0]] = line[1:-1]
                        else:
                            if float(line[11]) >= float(Blast_Dict[line[0]][10]):
                                Blast_Dict[line[0]] = line[1:-1]
                            elif float(line[11]) == float(Blast_Dict[line[0]][10]):
                                if randrange(0, 2) > 0:
                                    Blast_Dict[line[0]] = line[1:-1]
                            else:
                                pass
                    else:
                        pass

    # Convert dictionary to dataframe and export
    Blast_DF = pd.DataFrame.from_dict(Blast_Dict, orient='index')
    Blast_DF.to_csv(Output, sep='\t')

### ------------------------------- Main function ------------------------------
def main():
        parser = argparse.ArgumentParser(description='''Filters a Blast output either based on best bitscore (by default) of by both parameters provided and best bitscore'''
                                        'Global mandatory parameters: [Blast_File] [Output_File]\n'
                                        'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
        parser.add_argument('-b', '--blast', dest='Blast_File', action='store', required=True, help='Blast output of the FastA file search agains a DB')
        parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
        parser.add_argument('--id_perc', dest='ID_Perc', action='store', help='Minimum percentage identity for a match to be included. By default 40')
        parser.add_argument('--bistcore', dest='Bitscore', action='store', help='Minimum bitscore for a match to be included. By default 80')
        parser.add_argument('--evalue', dest='Evalue', action='store', help='Maximum Evalue for a match to be included. By default 0.1')
        parser.add_argument('--aln_percent', dest='Aln_Percent', action='store', help='If you have qlen and slen in your output, the minimum alignment the match must span to be included. By default not included')
        args = parser.parse_args()

        Blast_File = args.Blast_File
        Output_File = args.Output_File
        ID_Perc = args.ID_Perc
        Bitscore = args.Bitscore
        Evalue = args.Evalue
        Aln_Percent = args.Aln_Percent

        if ID_Perc == None:
            ID_Perc = 40
        if Bitscore == None:
            Bitscore = 80
        if Evalue == None:
            Evalue = 0.1

        Blast_Parser(Blast_File, Output_File, ID_Perc, Bitscore, Evalue, Aln_Percent)

if __name__ == "__main__":
    main()
