#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 26 March 2019

# Description: This script parses MAPLE's output and returns a table with each KEGG KO ID with its associated module and step in that module.
# MAPLE website: https://maple.jamstec.go.jp/maple/maple-2.3.1/index.html
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
import re
import argparse, sys
import pandas as pd

################################################################################

"""---2.0 Define Functions---"""

def Maple_Parser(Input_List, Annotation_File):
    KEGG_To_Module = pd.DataFrame(columns=['Gene','Module','Step','Total_Steps','Complex','FuncSet','Signature'])
    #KEGG_To_Module = {}
    for File in Input_List:
        with open(File) as Maple_out:
            Module = ""
            Steps = []
            Type = ""
            for line in Maple_out:
                line = line.strip()
                # Get module name and KO IDs involved in it.
                if re.match("^M.*\t", line):
                    Module = line.split(sep="\t")[0]
                    Steps = line.split(sep="\t")[2].split()
                    Type = line.split(sep="\t")[1]
                # Check if a KO ID is found in the annotation.
                elif re.match("K[0-9]{5}\t", line):
                    if int(line.split()[1]) > 0:
                        # Get KO ID
                        ID = line.split()[0]
                        # If KO not in the dataframe then create it with empty fields.
                        if ID not in KEGG_To_Module.index:
                            KEGG_To_Module.reindex(KEGG_To_Module.index.values.tolist()+[ID], fill_value="")
                            # If module is a pathway then append the module and the step in the pathway
                            if Type == "Pathway":
                                for index, value in enumerate(Steps):
                                    if ID in value:
                                        Index = index + 1
                                KEGG_To_Module.loc[ID, "Module"] = Module
                                KEGG_To_Module.loc[ID, "Step"] = str(Index)
                                KEGG_To_Module.loc[ID, "Total_Steps"] = str(len(Steps))
                            # If module is a complex then only append the module
                            elif Type == "Complex":
                                KEGG_To_Module.loc[ID, "Complex"] = Module
                            # If module is a funcset then only append the module
                            elif Type == "FuncSet":
                                KEGG_To_Module.loc[ID, "FuncSet"] = Module
                            elif Type == "Signature":
                                KEGG_To_Module.loc[ID, "Signature"] = Module
                            else:
                                print(File + " contains a type I dont recognize. Please check")
                                pass
                        elif ID in KEGG_To_Module.index:
                            if Type == "Pathway":
                                for index, value in enumerate(Steps):
                                    if ID in value:
                                        Index = index + 1
                                try:
                                    Updated_TotalSteps = KEGG_To_Module.loc[ID, "Total_Steps"] + ", " + str(len(Steps))
                                    Updated_module = KEGG_To_Module.loc[ID, "Module"] + ", " + Module
                                    Updated_Index = KEGG_To_Module.loc[ID, "Step"] + ", " + str(Index)
                                    KEGG_To_Module.loc[ID, "Module"] = Updated_module
                                    KEGG_To_Module.loc[ID, "Step"] = Updated_Index
                                    KEGG_To_Module.loc[ID, "Total_Steps"] = Updated_TotalSteps
                                except:
                                    KEGG_To_Module.loc[ID, "Module"] = Module
                                    KEGG_To_Module.loc[ID, "Step"] = str(Index)
                                    KEGG_To_Module.loc[ID, "Total_Steps"] = str(len(Steps))
                            # If module is a complex then only append the module
                            elif Type == "Complex":
                                try:
                                    KEGG_To_Module.loc[ID, "Complex"] = KEGG_To_Module.loc[ID, "Complex"] + ", " + Module
                                except:
                                    KEGG_To_Module.loc[ID, "Complex"] = Module
                            # If module is a funcset then only append the module
                            elif Type == "FuncSet":
                                try:
                                    KEGG_To_Module.loc[ID, "FuncSet"] = KEGG_To_Module.loc[ID, "FuncSet"] + ", " + Module
                                except:
                                    KEGG_To_Module.loc[ID, "FuncSet"] = Module
                            elif Type == "Signature":
                                try:
                                    KEGG_To_Module.loc[ID, "FuncSet"] = KEGG_To_Module.loc[ID, "FuncSet"] + ", " + Module
                                except:
                                    KEGG_To_Module.loc[ID, "Signature"] = Module
                            else:
                                print(File + " contains a type I dont recognize. Please check")
                                pass

    with open(Annotation_File) as Annotations:
        for line in Annotations:
            line = line.strip().split(sep="\t")
            if line[0] in KEGG_To_Module.index:
                KEGG_To_Module.loc[line[0], "Gene"] = line[1]

    return KEGG_To_Module

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''This script parses MAPLE's output and returns a table with each KEGG KO ID with its associated module and step in that module\n
                                                     The input should be any .matrix file from either the KAAS or the BLAST folders.'''
                                    'Global mandatory parameters: [Maple_Output] [Output Table]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--inputMaple", dest='Input_Mapple', action='store', required=True, nargs='+', help="Input Maple Annotation")
    parser.add_argument('-a', '--annotations', dest='Annotation_File', action='store', required=True, help='Mapping file from KO IDs to gene, downloaded from https://www.genome.jp/kegg-bin/get_htext#B1')
    parser.add_argument('-o', '--outputTable', dest='Output_Table', action='store', required=True, help='Output table in the form KEGG_ID Module Step')
    args = parser.parse_args()

    Input_Mapple = args.Input_Mapple
    Annotation_File = args.Annotation_File
    Output_Table = args.Output_Table


    # Run parser
    Parsed_Mapple = Maple_Parser(Input_Mapple, Annotation_File)
    Parsed_Mapple.to_csv(Output_Table, sep="\t")

if __name__ == "__main__":
    main()
