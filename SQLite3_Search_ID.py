#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  2.0
# Date:		 August 12 2019

# Description: This searches IDs provided in a given database.
# Its use is mainly to search for RefSeq and Uniprot IDs faster than parsing files.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""


def Search_Table(SQL_database, Input_List, Database, Output, Columns = None, Standalone = True):
    import pandas as pd
    import sqlite3
    
    conn = sqlite3.connect(SQL_database)
    cur = conn.cursor()

    # Create empty list to store annotations.
    Annotation_List = []

    # Search the DB and append results to Annotation_List
    print("\t".join(Columns))
    #with open(Output) as Final_Table:

    if Standalone:
        # Set column names.
        if Database == "Swissprot":
            Col_Names = ['ID' , 'Accession', 'Name', 'KO_Uniprot', 'Organism', 'Taxonomy', 'Function', 'Compartment', 'Process']
        elif Database == "Trembl":
            Col_Names = ['ID' , 'Accession', 'Name', 'KO_Uniprot', 'Organism', 'Taxonomy', 'Function', 'Compartment', 'Process']
        elif Database == "RefSeq":
            Col_Names = ['ID', 'Gene_Name', 'Taxonomy', 'Note']

        for ID in (Input_List):
            cur.execute("SELECT * FROM " + Database + " WHERE ID=?", (ID,))
            rows = cur.fetchall()
            try:
                Annotation_List.append(rows[0])
            except:
                pass
        Annotation_DF = pd.DataFrame(Annotation_List, columns=Col_Names)
        Annotation_DF.set_index('ID', inplace=True)
        Annotation_DF.to_csv(Output, sep="\t")

    else:
        for ID in (Input_List):
            cur.execute("SELECT * FROM " + Database + " WHERE ID=?", (ID,))
            rows = cur.fetchall()
            try:
                Annotation_List.append(rows[0])
            except:
                pass
        return Annotation_List

    cur.close()
    conn.close()


################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script builds a sqlite database from a tab separated table.\n'''
            '''By default it assumes the first line of the input table has headers, if not\n'''
            '''you must provide a list of headers as -h header1,header2,header3...\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input Table] -d [Database Name] -t [Table Name] -n [Index Name] --header [Header List]\n'''
            '''Global mandatory parameters: -i [Input Table] -d [Database Name] -t [Table Name] -n [Index Name]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Input_Table', action='store', required=True, help='Input tab-delimited table to parse, by default assumes headers are present')
    parser.add_argument('-o', '--output', dest='Output_Table', action='store', required=True, help='Output tab-delimited table to store annotations')
    parser.add_argument('-s', '--sql_database', dest='SQL_database', action='store', required=True, help='SQL database where the annotations are stored')
    parser.add_argument('--gene_col', dest='Gene_Col', action='store', required=False, type=int,
                        help='Column with gene IDs in the input file, by default None, i.e. assumes only IDs of hits are given.')
    parser.add_argument('--id_col', dest='ID_Col', action='store', required=True, type=int, 
                        help='Column with IDs of database hits in the input file')
    parser.add_argument('--database', dest='Database', action='store', required=False, default='Swissprot',
                        help='Comma-separated names of databases to search, can include either "Swissprot", "Trembl", or "RefSeq". By default "Swissprot"')
    args = parser.parse_args()

    Input_Table = args.Input_Table
    Output_Table = args.Output_Table
    SQL_database = args.SQL_database
    Gene_Col = args.Gene_Col
    ID_Col = args.ID_Col
    Database = args.Database

    # Set column names.
    if Database == "Swissprot":
        Col_Names = ['ID' , 'Accession', 'Name', 'KO_Uniprot', 'Organism', 'Taxonomy', 'Function', 'Compartment', 'Process']
    elif Database == "Trembl":
        Col_Names = ['ID' , 'Accession', 'Name', 'KO_Uniprot', 'Organism', 'Taxonomy', 'Function', 'Compartment', 'Process']
    elif Database == "RefSeq":
        Col_Names = ['ID', 'Gene_Name', 'Taxonomy', 'Note']

    # Create table with or without gene ids
    if Gene_Col != None:
        Col_Names.insert(1,'Gene_ID')
        Annotation_Dict = {}
        Input_List = []
        with open(Input_Table) as ID_File:
            for line in ID_File:
                line = line.strip().split()
                Gene = line[Gene_Col - 1]
                if Database == "Swissprot" or Database == "Trembl":
                    Hit = line[ID_Col - 1].split("|")[2]
                else:
                    Hit = line[ID_Col - 1]
                Annotation_Dict[Hit] = [Gene]
                Input_List.append(Hit)
        # Get list of annotations
        Annotation_List = Search_Table(SQL_database, Input_List, Database, Standalone=False)
        for record in Annotation_List:
            if len(Annotation_Dict[record[0]]) > 2:
                pass
            else:
                Annotation_Dict[record[0]].extend(list(record[1:]))

        Annotation_DF = pd.DataFrame.from_dict(Annotation_Dict, orient='index')
        Annotation_DF.reset_index(inplace=True)
        Annotation_DF.columns = Col_Names
        Annotation_DF.set_index('Gene_ID', inplace=True)
        Annotation_DF.to_csv(Output_Table, sep="\t")
    else:
        Input_List = []
        with open(Input_Table) as ID_File:
            for line in ID_File:
                line = line.strip().split()
                if Database == "Swissprot" or Database == "Trembl":
                    Hit = line[ID_Col - 1].split("|")[2]
                else:
                    Hit = line[ID_Col - 1]
                Input_List.append(Hit)
        Annotation_List = Search_Table(SQL_database, Input_List, Database, Output=Output_Table, Standalone=True)


if __name__ == "__main__":
    main()