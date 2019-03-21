#!/usr/bin/env python

################################################################################
"""---1.0 Import Modules---"""

import re
import pandas as pd
import argparse, sys

################################################################################
"""---2.0 Define Functions---"""

def ParseClusterFile(ClusterFile):
    # Create empty Cluster variable, Empty dictionary and empty DataFrame.
    regex = re.compile("_[0-9]+_[0-9]+\.\.\.")
    Cluster = None
    Clustr_Gene = {}
    PangenomeMatrix = pd.DataFrame()
    ClusterRep = pd.DataFrame(columns=['Cluster','Representative_Gene'])
    # Open file and iterate through lines
    with open(ClusterFile) as file:
        for line in file:
            line = line.strip().replace(">", "").split()
            if len(line) == 2: # Looking for "Cluster #"
                # Add Cluster ID to dictionary and as index in the Dataframe.
                Cluster = line[0] + " " + line[1]
                Clustr_Gene[Cluster] = []
                PangenomeMatrix = PangenomeMatrix.reindex(PangenomeMatrix.index.values.tolist()+[Cluster], fill_value=0)
            else:
                # Get Gene and Genome names and append the gene to the dictionary.
                Gene = line[2].replace("...", "")
                Genome = re.sub(regex, "", line[2])
                Clustr_Gene[Cluster].append(Gene)
                if line[3] == "*":
                    ClusterRep = ClusterRep.append({'Cluster' : Cluster , 'Representative_Gene' : Gene} , ignore_index=True)
                # Insert Genome in dataframe and insert gene at position "Cluster #", "Genome".
                if Genome not in PangenomeMatrix:
                    PangenomeMatrix[Genome] = 0
                PangenomeMatrix.at[Cluster, Genome] = 1

    # Return both dataframes
    return PangenomeMatrix, ClusterRep

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
    parser.add_argument('-r', '--representative', dest='Representative_File', action='store', help='Output file for the representatives table, if none, "Representatives.tab"', default="Representatives.tab")
    args = parser.parse_args()

    Cluster_File = args.Cluster_File
    Output_File = args.Output_File
    Representative_File = args.Representative_File

    # Return dataframes
    PangenomeMatrix, ClusterRep = ParseClusterFile(Cluster_File)

    # Export both dataframes to tab-separated tables.
    PangenomeMatrix.to_csv(Output_File, sep='\t')
    ClusterRep.to_csv(Representative_File, sep='\t', index=False)

if __name__ == "__main__":
    main()
