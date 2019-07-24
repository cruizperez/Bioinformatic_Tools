#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Institution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 24 July 2019

# Description: This script parses a Uniprot.dat file and outputs a table with
# the ID, Accession, Gene Name, Organism, Taxonomy, KEGG ID, Function, Compartment, and Process.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""

from collections import defaultdict
import pandas as pd

################################################################################
"""---2.0 Define Functions---"""

def Parse_Uniprot(Uniprot_Dat):
    Uniprot_Dictionary = defaultdict(list)
    with open(Uniprot_Dat) as Uniprot:
        for line in Uniprot:
            if line.startswith("ID"):
                ID = line.split()[1]
                Uniprot_Dictionary[ID] = [[] for k in range(8)]
            elif line.startswith("AC"):
                Uniprot_Dictionary[ID][0].append(line.split()[1].replace(";", ""))
            elif "RecName" in line:
                Name = line.split("Full=")[1].replace(";", "")
                Name = Name.split("{")[0].strip()
                Uniprot_Dictionary[ID][1].append(Name)
            elif line.startswith("OS"):
                Uniprot_Dictionary[ID][2].append(line.split("OS")[1].strip())
            elif line.startswith("OC"):
                Uniprot_Dictionary[ID][3].append(line.split("OC")[1].strip())
            elif line.startswith("DR"):
                if "KO;" in line:
                    Uniprot_Dictionary[ID][4].append(line.split()[2].replace(";", ""))
                elif "F:" in line:
                    Uniprot_Dictionary[ID][5].append(line.split("F:")[1].split(";")[0])
                elif "C:" in line:
                    Uniprot_Dictionary[ID][6].append(line.split("C:")[1].split(";")[0])
                elif "P:" in line:
                    Uniprot_Dictionary[ID][7].append(line.split("P:")[1].split(";")[0])
    return Uniprot_Dictionary

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''This script parses a Uniprot.dat file and outputs a table with\n'''
                                                    '''the ID, Accession, Gene Name, Organism, Taxonomy, KEGG ID, Function, Compartment, and Process.\n'''
                                    '''\nGlobal mandatory parameters: [Input Uniprot.dat File]\n'''
                                    '''\nOptional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Uniprot_File', action='store', required=True, help='Uniprot.dat file to parse')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=False, help='Output table, if none "Uniprot_Accessions.tab".', default="Uniprot_Accessions.tab")
    args = parser.parse_args()

    Uniprot_File = args.Uniprot_File
    Output_File = args.Output_File

    Dictionary = Parse_Uniprot(Uniprot_File)
    Uniprot_DF = pd.DataFrame.from_dict(Dictionary,orient='index', dtype="str")
    Uniprot_DF.columns = ["Accesion", "Gene", "Organism", "Taxonomy", "KEGG", "Function", "Compartment", "Process"]
    for column in Uniprot_DF:
        Uniprot_DF[column] = Uniprot_DF[column].str[0]
        Uniprot_DF.to_csv(Output_File, sep="\t")

if __name__ == "__main__":
    main()
