#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Institution:   Georgia Institute of Technology
# Version:	  1.1
# Date:		 29 July 2019

# Description: This script parses a Uniprot.dat file and outputs a table with
# the ID, Accession, Gene Name, KO NUmber, Organism, Taxonomy, KEGG ID, Function, Compartment, and Process.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Parse_Uniprot(Uniprot_Dat, Output, Header = False):
    Output = open(Output, 'w')
    if Header == True:
        Output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ID", "Accession", "Gene_Name", "KO Number", "Organism", "Taxonomy", "Function", "Compartment", "Process"))
    with open(Uniprot_Dat) as Uniprot:
        ID = ""
        Accession = ""
        Name = ""
        KO = ""
        Organism = ""
        Taxonomy = ""
        Function = ""
        Compartment = ""
        Process = ""
        for line in Uniprot:
            if line.startswith("ID", 0):
                ID = line.split()[1]
            elif "AC  " in line:
                Accession = line.split()[1]
            elif "RecName" in line:
                Name = line.split("Full=")[1]
                Name = Name.split("{")[0].strip()
            elif "OS  " in line:
                Organism = ' '.join([Organism, line.split("OS")[1].strip()])
            elif "OC  " in line:
                Taxonomy = ' '.join([Taxonomy, line.split("OC")[1].strip()])
            elif "DR   KO;" in line:
                KO = line.split()[2].replace(";", "")
            elif "DR   GO;" in line:
                if "; F:" in line:
                    Function = ''.join([Function, line.split("GO;")[1].strip(), " -- "])
                elif "; C:" in line:
                    Compartment = ''.join([Compartment, line.split("GO;")[1].strip(), " -- "])
                elif "; P:" in line:
                    Process = ''.join([Process, line.split("GO;")[1].strip(), " -- "])
            elif "//\n" in line:
                Output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ID, Accession, Name, KO, Organism, Taxonomy, Function, Compartment, Process))
                ID = ""
                Accession = ""
                Name = ""
                KO = ""
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

    Parse_Uniprot(Uniprot_File, Output_File, Headers)

if __name__ == "__main__":
    main()
