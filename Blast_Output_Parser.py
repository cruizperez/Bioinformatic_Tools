#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import sys, argparse, os
import pandas as pd
from random import randrange

"""----------------------------- 1.0 Define Functions -----------------------------"""
### ------------------------------Match filter--------------------------------------
def HitConfidence(line, id, bitscore, evalue, Aln_Percent = None, Shorter = None, Query = None, Subject = None):
    if Aln_Percent == None:
        if float(line[2]) >= float(id) and float(line[11]) >= float(bitscore) and float(line[10]) <= float(evalue):
            High_Quality_Match = True
        else:
            High_Quality_Match = False
    else:
        if Shorter == True:
            if (int(line[3])*100/min(int(line[13]),int(line[12])) >= float(Aln_Percent)) and float(line[11]) >= float(bitscore) and float(line[2]) >= float(id) and float(line[10]) <= float(evalue):
                High_Quality_Match = True
            else:
                High_Quality_Match = False
        elif Query == True:
            if (int(line[3])*100/int(line[12])) >= float(Aln_Percent) and float(line[11]) >= float(bitscore) and float(line[2]) >= float(id) and float(line[10]) <= float(evalue):
                High_Quality_Match = True
            else:
                High_Quality_Match = False
        elif Subject == True:
            if (int(line[3])*100/int(line[13])) >= float(Aln_Percent) and float(line[11]) >= float(bitscore) and float(line[2]) >= float(id) and float(line[10]) <= float(evalue):
                High_Quality_Match = True
            else:
                High_Quality_Match = False
        else:
            if (int(line[3])*100/max(int(line[13]),int(line[12])) >= float(Aln_Percent)) and float(line[11]) >= float(bitscore) and float(line[2]) >= float(id) and float(line[10]) <= float(evalue):
                High_Quality_Match = True
            else:
                High_Quality_Match = False
    return(High_Quality_Match)

### ----------------------------- Blast Parser ----------------------------------------
def Blast_Parser(BlastFile, Output, id, bitscore, evalue, Aln_Percent = None, Shorter = None, Query = None, Subject = None):
    Blast_Dict = {}
    print("Reading " + BlastFile + " Blast Output")
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
        print("Only retrieving best match per sequence using other default values.")
        if id == None:
            id = 30
        if bitscore == None:
            bitscore = 50
        if evalue == None:
            evalue = 10
        # Give only best match based on bitscore.
        with open(BlastFile) as BlastFile_Input:
                for line in BlastFile_Input:
                    line = line.strip().split("\t")
                    Good = HitConfidence(line, id, bitscore, evalue, Aln_Percent, Shorter, Query, Subject)
                    if Good == True:
                        if line[0] not in Blast_Dict:
                            Blast_Dict[line[0]] = line[1:]
                        else:
                            if float(line[11]) >= float(Blast_Dict[line[0]][10]):
                                Blast_Dict[line[0]] = line[1:]
                            elif float(line[11]) == float(Blast_Dict[line[0]][10]):
                                    if randrange(0, 2) > 0:
                                        Blast_Dict[line[0]] = line[1:]
                            else:
                                pass
                    else:
                        pass
    else:
        print("Performing match filtering and retrieving best match based on the parameters you provided (if only some, the others will take their default values).")
        if id == None:
            id = 30
        if bitscore == None:
            bitscore = 50
        if evalue == None:
            evalue = 10
        # Do match filtering based on parameters provided and retrieve best match based on best bitscore.
        with open(BlastFile) as BlastFile_Input:
            if long == True:
                for line in BlastFile_Input:
                    line = line.strip().split("\t")
                    Good = HitConfidence(line, id, bitscore, evalue, Aln_Percent, Shorter, Query, Subject)
                    if Good == True:
                        if line[0] not in Blast_Dict:
                            Blast_Dict[line[0]] = line[1:]
                        else:
                            if float(line[11]) >= float(Blast_Dict[line[0]][10]):
                                Blast_Dict[line[0]] = line[1:]
                            elif float(line[11]) == float(Blast_Dict[line[0]][10]):
                                if randrange(0, 2) > 0:
                                    Blast_Dict[line[0]] = line[1:]
                            else:
                                pass
                    else:
                        pass
            else:
                for line in BlastFile_Input:
                    line = line.strip().split("\t")
                    Good = HitConfidence(line, id, bitscore, evalue)
                    if Good == True:
                        if line[0] not in Blast_Dict:
                            Blast_Dict[line[0]] = line[1:]
                        else:
                            if float(line[11]) >= float(Blast_Dict[line[0]][10]):
                                Blast_Dict[line[0]] = line[1:]
                            elif float(line[11]) == float(Blast_Dict[line[0]][10]):
                                if randrange(0, 2) > 0:
                                    Blast_Dict[line[0]] = line[1:]
                            else:
                                pass
                    else:
                        pass

    # Convert dictionary to dataframe and export
    Blast_DF = pd.DataFrame.from_dict(Blast_Dict, orient='index')
    Blast_DF.to_csv(Output, sep='\t', header= False)

### ------------------------------- Main function ------------------------------
def main():
        parser = argparse.ArgumentParser(description='''Filters a Blast output either based on best bitscore (by default) of by both parameters provided and best bitscore'''
                                        'Global mandatory parameters: [Blast_File] [Output_File]\n'
                                        'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
        parser.add_argument('-b', '--blast', dest='Blast_File', action='store', required=True, help='Blast output of the FastA file search agains a DB')
        parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
        parser.add_argument('--id_perc', dest='ID_Perc', action='store', help='Minimum percentage identity for a match to be included. By default 30')
        parser.add_argument('--bistcore', dest='Bitscore', action='store', help='Minimum bitscore for a match to be included. By default 50')
        parser.add_argument('--evalue', dest='Evalue', action='store', help='Maximum Evalue for a match to be included. By default 10')
        parser.add_argument('--aln_percent', dest='Aln_Percent', action='store', help='If you have qlen and slen in your output, the minimum alignment the match must span to be included. By default not included')
        parser.add_argument('--shorter', action='store_true', help='Calculates the alignment percentage on the shorter sequence, by default false, i.e. calculated on the longer')
        parser.add_argument('--query', action='store_true', help='Calculates the alignment percentage on the query sequence, by default false, i.e. calculated on the longer')
        parser.add_argument('--subject', action='store_true', help='Calculates the alignment percentage on the subject sequence, by default false, i.e. calculated on the longer')
        args = parser.parse_args()

        Blast_File = args.Blast_File
        Output_File = args.Output_File
        ID_Perc = args.ID_Perc
        Bitscore = args.Bitscore
        Evalue = args.Evalue
        Aln_Percent = args.Aln_Percent
        Shorter = args.shorter
        Query = args.query
        Subject = args.subject

        Blast_Parser(Blast_File, Output_File, ID_Perc, Bitscore, Evalue, Aln_Percent, Shorter, Query, Subject)

if __name__ == "__main__":
    main()
