#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 22 March 2019

# Description: This script parses Blast 2.2.31+ SAM output )-outfmt 15), selecting the best hit
and removing the duplications in the reference IDs
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""
from random import randrange
import argparse, sys

################################################################################
"""---2.0 Define Functions---"""

def Blast_SAM_Parser(SAM_File, Output_File):
    Output_FH = open(Output_File, "w")
    with open(SAM_File, "r") as SAM_FH:
        Ref_ID = []
        Read_ID = {}
        for line in SAM_FH:
            line = line.strip()
            if line.startswith('@SQ'):
                Ref = line.split()[1]
                Ref = Ref.split("|")[1]
                if Ref not in Ref_ID:
                    Ref_ID.append(Ref)
                    Output_FH.write(line + "\n")
                else:
                    pass
            elif line.startswith('@'):
                pass
                Output_FH.write(line + "\n")
            else:
                Read = line.split()[0]
                Bitscore = float(line.split("BS:f:")[1])
                if Read not in Read_ID:
                    Read_ID[Read] = [Bitscore, line]
                else:
                    if Bitscore > Read_ID[Read][0]:
                        Read_ID[Read] = [Bitscore, line]
                    elif Bitscore == Read_ID[Read][0]:
                        if randrange(0,2) > 0:
                            Read_ID[Read] = [Bitscore, line]
                        else:
                            pass
                    else:
                        pass

        for key in Read_ID:
            Output_FH.write(Read_ID[key][1] + "\n")


################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''This script parses Blast 2.2.31+ SAM output (-outfmt 15), selecting the best hit
                                                    \nand removing the duplications in the reference IDs\n'''
                                    'Global mandatory parameters: [SAM_File] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-s', '--sam', dest='SAM_File', action='store', required=True, help='SAM file to parse')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output SAM file best hits')
    args = parser.parse_args()

    SAM_File = args.SAM_File
    Output_File = args.Output_File

    Blast_SAM_Parser(SAM_File, Output_File)

if __name__ == "__main__":
    main()
