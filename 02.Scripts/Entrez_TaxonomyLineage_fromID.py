#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 August 4 2019

# Description: This script takes an NCBI ID of a genome and returns the taxonomy ID.
If you provide a list instead it will give you a tab-separated list with [Genome ID] [Taxonomy ID] [Lineage]
########################################################################
"""
################################################################################
"""---1.0 Import Modules---"""
from Bio import Entrez

################################################################################
"""---2.0 Define Functions---"""
def Post_IDs(id_list, Database, email, api_key=None):
    if Database == "biosample":
        import re
        for index, item in enumerate(id_list):
            id_list[index] = re.sub('[A-Z]+', '', item)
    Entrez.api_key = api_key
    Entrez.email = email
    search_results = Entrez.read(Entrez.epost(Database, id=",".join(id_list)))
    print(search_results)
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    History = (webenv, query_key)

    return History

#----------------------------------------

def Find_Taxonomy_FromID(ID_List, Database, History, email, API= None):
    Batch_Size = 2000
    Count = len(ID_List)
    Entrez.api_key = API
    Entrez.email = email
    Taxonomy_Dictionary = {}
    #Download in batches
    for start in range(0, Count, Batch_Size):
        end = min(Count, start+Batch_Size)
        print("Going to download record {} to {}".format(start+1, end))
        fetch_handle = Entrez.esummary(db=Database,
                                    rettype="native", retmode="text",
                                    retstart=start, retmax=Batch_Size,
                                    webenv=History[0], query_key=History[1])
        Data_Dicts = Entrez.read(fetch_handle)
        for Dict in Data_Dicts:
            TaxID = Dict['TaxId']
            Title = Dict['Title']
            AccessionVersion = Dict['AccessionVersion']
            Taxonomy_Dictionary[AccessionVersion] = [str(TaxID), Title]
        fetch_handle.close()

    return Taxonomy_Dictionary

# -------------------------------------------

def Get_Lineage_FromTaxID(TaxID_List, History, email, API = None, Other = None):
    Batch_Size = 500
    Count = len(TaxID_List)
    Entrez.api_key = API
    Entrez.email = email
    Original_Dict = Other
    for start in range(0, Count, Batch_Size):
        end = min(Count, start+Batch_Size)
        print("Going to download record {} to {}".format(start+1, end))
        fetch_handle = Entrez.efetch(db="taxonomy", retmode="xml",
                                     retstart=start, retmax=Batch_Size,
                                     webenv=History[0], query_key=History[1])
        Data_Dicts = Entrez.read(fetch_handle)
        if Original_Dict != None:
            for Dict in Data_Dicts:
                Tax_ID = Dict['TaxId']
                Lineage = Dict['Lineage']
                for key, value in Original_Dict.items():
                    if value[0] == Tax_ID and len(value) < 3:
                        Original_Dict[key].append(Lineage)
            return Original_Dict
        else:
            Lineage_Dictionary = {}
            for Dict in Data_Dicts:
                Tax_ID = Dict['TaxId']
                Lineage = Dict['Lineage']
                Lineage_Dictionary[Tax_ID] = Lineage
            return Lineage_Dictionary
        fetch_handle.close()

#--------------------------------------

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script takes an NCBI ID of a genome and returns the taxonomy ID.\n'''
            '''If you provide a list instead it will give you a tab-separated list\n'''
            '''with [Genome ID] [Taxonomy ID] [Lineage]\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Genome ID] or -l [ID list file] -e [Email to tell NCBI who you are. IMPORTANT] -o [Output_Table (For list)] -a [API Key obtained from NCBI]\n'''
            '''Global mandatory parameters: -e [Email]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-d', '--dbfrom', dest='DB_From', action='store', required=False, default="nuccore", help='Database to look from, by default "nuccore"')
    parser.add_argument('-t', '--targetdb', dest='DB_Target', action='store', required=False, default="taxonomy", help='Database to look from, by default "nuccore"')
    parser.add_argument('--item', dest='Item', action='store', required=False, default="id", help='Item to recover, e.g., id (default) or lineage.')
    parser.add_argument('-id_list', '--id_list', dest='id_list',
        action='store', required=False, nargs="+",
        help='Comma-separated list of IDs, e.g., ID1,ID2,ID3...')
    parser.add_argument('-f', '--fileID', dest='ID_File', action='store', required=False, help='File with IDs to search, one per line')
    parser.add_argument('-o', '--output', dest='Output_Table', action='store', required=True, help='Table to which the taxonomy and lineages will be printed.')
    parser.add_argument('-e', '--email', dest='Email', action='store', required=True, help='Email to tell NCBI who you are.')
    parser.add_argument('-a', '--api', dest='API', action='store', required=False, help='API key for large queries')
    parser.add_argument('--header', dest='Header', action='store_true', required=False, help='Add header to output, False by default')
    args = parser.parse_args()

    DB_From = args.DB_From.lower()
    DB_Target = args.DB_Target.lower()
    Item = args.Item.lower()
    id_list = args.id_list
    ID_File = args.ID_File
    Output_Table = args.Output_Table
    Email = args.Email.lower()
    API = args.API.lower()
    Header = args.Header

    Input_ID_List = []
    if ID_File == None and id_list == None:
        print("No IDs provided. Please give a list of IDs with -l or a file containing the IDs one per line with -f.")
    elif ID_File != None:
        with open(ID_File) as Input_File:
            for line in Input_File:
                line = line.strip().split()[0]
                Input_ID_List.append(line)
    elif id_list != None:
        Input_ID_List = id_list
    print(Input_ID_List)
    print("Performing the search of {} records in {} database and returning {} from {} database".format(len(Input_ID_List),
    DB_From, Item, DB_Target))

    if DB_From == "nuccore" and DB_Target == "taxonomy":
        # First history search in nuccore
        Nuccore_History = Post_IDs(Input_ID_List, DB_From, Email, API)
        # Retrieve the IDs from the target database
        NucID_TaxID_Dict = Find_Taxonomy_FromID(Input_ID_List, DB_From, Nuccore_History, Email, API)
        if Item == "id":
            if Header == False:
                ID_Dataframe = pd.DataFrame.from_dict(NucID_TaxID_Dict, orient='index')
            else:
                ID_Dataframe = pd.DataFrame.from_dict(NucID_TaxID_Dict, orient='index', columns=["TaxID", "Element"])
                ID_Dataframe.index.name = 'ID'
            ID_Dataframe.to_csv(Output_Table, sep="\t", index=True)
        elif Item == "lineage":
            Input_Tax_IDs = []
            for key, value in NucID_TaxID_Dict.items():
                Input_Tax_IDs.append(value[0])
            # Second search history
            Taxonomy_History = Post_IDs(Input_Tax_IDs, DB_Target, Email, API)
            # Search for lineages based on TaxId
            Lineage_Dictionary = Get_Lineage_FromTaxID(Input_Tax_IDs, Taxonomy_History, Email, API, NucID_TaxID_Dict)
            if Header == False:
                Lineage_Dataframe = pd.DataFrame.from_dict(Lineage_Dictionary, orient='index')
            else:
                Lineage_Dataframe = pd.DataFrame.from_dict(Lineage_Dictionary, orient='index', columns=["TaxID", "Element", "Lineage"])
                Lineage_Dataframe.index.name = 'ID'
            Lineage_Dataframe.to_csv(Output_Table, sep="\t", index=True)


if __name__ == "__main__":
    main()
