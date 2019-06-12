#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 05 June 2019

# Description: This script parses a CD-HIT or a MeShClust Cluster File and returns
the longest sequence per cluster.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""
from pandas import DataFrame
import argparse, sys


################################################################################
"""---2.0 Define Functions---"""

def LongestExtract(Input_ClusterFile):
    Cluster_Dict = {}
    Len2Add = 0
    with open(Input_ClusterFile) as ClusterFile:
        for line in ClusterFile:
            line = line.strip()
            if line.startswith(">Cluster"):
                ClusterID = line.split('>')[1]
                Len2Add = 0
            else:
                line = line.split('\t')[1].split(' ')
                ContigLen = int(line[0].replace('nt,', ''))
                if ContigLen > Len2Add:
                    ContigName = line[1].replace('>', '').replace('...','')
                    Cluster_Dict[ClusterID] = ContigName
                else:
                    continue
    return Cluster_Dict

################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''This script parses a CD-HIT or a MeShClust Cluster File
                                                    \nand returns the longest sequence per cluster.'''
                                    'Global mandatory parameters: [Input_Cluster_File] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Input_ClusterFile', action='store', required=True, help='Cluster file to parse')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output Table')
    args = parser.parse_args()

    Input_ClusterFile = args.Input_ClusterFile
    Output_File = args.Output_File

    Cluster_Dictionary = LongestExtract(Input_ClusterFile)
    Cluster_DF = DataFrame.from_dict(Cluster_Dictionary, orient='index')
    Cluster_DF.to_csv(Output_File, sep='\t', header= False)

if __name__ == "__main__":
    main()
