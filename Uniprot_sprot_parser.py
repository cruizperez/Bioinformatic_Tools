#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Institution:   Georgia Institute of Technology
# Version:	  1.1
# Date:		 29 July 2019

# Description: This script parses a Uniprot.dat file and outputs a table with
# the ID, Accession, Gene Name, Organism, Taxonomy, KEGG ID, Function, Compartment, and Process.
# To use in a faster way use gnu parallel as follows:
# cat InputFile.dat | parallel --jobs [#] --pipe --recend '//' cat \> File_{#}\;
# /mnt/c/Users/Cruiz/Documents/GitHub/Misc_Tools/Uniprot_sprot_parser.py -i File_{#} -o Test.parsed\; rm File_{#}\;
# This will split the file into records, execute the script and remove the splitted file. Remember to give it a
# single output so everything will be in the same file.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""

from collections import defaultdict
import pandas as pd

################################################################################
"""---2.0 Define Functions---"""

def Parse_Uniprot(Uniprot_Dat, Output, Header = False):
    Uniprot_Dictionary = defaultdict(list)
    Output = open(Output, 'w')
    if Header == True:
        Output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ID", "Accession", "Gene_Name", "Organism", "Taxonomy", "Function", "Compartment", "Process"))
    with open(Uniprot_Dat) as Uniprot:
        ID = ""
        Accession = ""
        Name = ""
        Organism = ""
        Taxonomy = ""
        Function = ""
        Compartment = ""
        Process = ""
        for line in Uniprot:
            if "ID  " in line:
                ID = line.split()[1]
                Uniprot_Dictionary[ID] = [[] for k in range(8)]
            elif "AC  " in line:
                Accession = line.split()[1]
                Uniprot_Dictionary[ID][0].append(line.split()[1].replace(";", ""))
            elif "RecName" in line:
                Name = line.split("Full=")[1]
                Name = Name.split("{")[0].strip()
                #Uniprot_Dictionary[ID][1].append(Name)
            elif "OS  " in line:
                Organism = ' '.join([Organism, line.split("OS")[1].strip()])
            elif "OC  " in line:
                Taxonomy = ' '.join([Taxonomy, line.split("OC")[1].strip()])
            elif "DR  " in line:
                if "KO;" in line:
                    KO = line.split()[2]
                elif "; F:" in line:
                    Function = ''.join([Function, line.split("GO;")[1].strip(), " -- "])
                elif "; C:" in line:
                    Compartment = ''.join([Compartment, line.split("GO;")[1].strip(), " -- "])
                elif "; P:" in line:
                    Process = ''.join([Process, line.split("GO;")[1].strip(), " -- "])
            elif "//\n" in line:
                Output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ID, Accession, Name, Organism, Taxonomy, Function, Compartment, Process))
                ID = ""
                Accession = ""
                Name = ""
                Organism = ""
                Taxonomy = ""
                Function = ""
                Compartment = ""
                Process = ""

    Output.close()

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''This script parses a Uniprot.dat file and outputs a table with\n'''
                                                    '''the ID, Accession, Gene Name, Organism, Taxonomy, KEGG ID, Function, Compartment, and Process.\n
                                                    For faster usage in alrge files use gnu parallel (read script file to see how)\n'''
                                    '''\nGlobal mandatory parameters: [Input Uniprot.dat File]\n'''
                                    '''\nOptional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Uniprot_File', action='store', required=True, help='Uniprot.dat file to parse')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=False, help='Output table')
    parser.add_argument('--header', dest='Headers', action='store_true', required=False, help='Output table should contain headers, False by default')
    args = parser.parse_args()

    Uniprot_File = args.Uniprot_File
    Output_File = args.Output_File
    Headers = args.Headers

    # Create empty dataframe with colnames.
    #with open(Output_File, 'w') as OutFile:
    #    OutFile.write("ID\tAccesion\tGene\tOrganism\tTaxonomy\tKEGG\tFunction\tCompartment\tProcess\n")
    # Get dictionary and convert to df
    Dictionary = Parse_Uniprot(Uniprot_File, Output_File)

if __name__ == "__main__":
    main()
