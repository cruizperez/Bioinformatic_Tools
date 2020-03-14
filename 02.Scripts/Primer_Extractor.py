#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 24 March 2019

# Description: This script splits a FastA file into a fiven number of files.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""

import re
import argparse, sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

################################################################################
"""---2.0 Define Functions---"""

def Primer_Extractor(FastA_Input, Output_File, PrimerF_List = None, PrimerR_List = None):
    # Create empty Cluster variable, Empty dictionary and empty DataFrame.
    Output = open(Output_File, 'w')
    with open(FastA_Input) as FastA_File:
        for title, seq in SimpleFastaParser(FastA_File):
            for PrimerF in PrimerF_List:
                for PrimerR in PrimerR_List:
                    #print(PrimerF)
                    #PrimerF = PrimerF.replace("T", "U")
                    #PrimerR = PrimerR.replace("T", "U")
                    try:
                        Start = re.search(PrimerF, seq).start()
                    except:
                        Start = None
                    try:
                        End = re.search(PrimerR, seq).start()
                    except:
                        End = None
                    if End != None:
                        Output.write(">%s\n%s\n" % (title, seq[0:-1]))


################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Parses a CD-HIT-type cluster file and returns a matrix for cluster presence and
                                                    a table with the representative gene per cluster'''
                                    'Global mandatory parameters: [Cluster File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Cluster_File', action='store', required=True, help='Clustr file to parse')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output file for the matrix counts, if none "Matrix_Counts.tab".', default="Matrix_Counts.tab")
    parser.add_argument('-f', '--forward', dest='Forward_Primer', nargs = '+', help='Output file for the representatives table, if none, "Representatives.tab"', default="Representatives.tab")
    parser.add_argument('-r', '--reverse', dest='Reverse_Primer', nargs = '+', help='Output file for the representatives table, if none, "Representatives.tab"', default="Representatives.tab")
    args = parser.parse_args()

    Cluster_File = args.Cluster_File
    Output_File = args.Output_File
    Forward_Primer = args.Forward_Primer
    Reverse_Primer = args.Reverse_Primer

    Primer_Extractor(Cluster_File, Output_File, Forward_Primer, Reverse_Primer)
    # Export both dataframes to tab-separated tables.
    #PangenomeMatrix.to_csv(Output_File, sep='\t')
    #ClusterRep.to_csv(Representative_File, sep='\t', index=False)

if __name__ == "__main__":
    main()
