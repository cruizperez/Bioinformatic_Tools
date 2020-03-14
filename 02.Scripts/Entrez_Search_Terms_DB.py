#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  2.0
# Date:		 November 6 2019

# Description: This script takes an NCBI ID of a genome and returns the taxonomy ID.
If you provide a list instead it will give you a tab-separated list with [Genome ID] [Taxonomy ID] [Lineage]
########################################################################
"""
################################################################################
"""---1.0 Define Functions---"""
def Post_IDs(List, Database, email, API= None):
    from Bio import Entrez
    if Database == "biosample":
        import re
        for index, item in enumerate(List):
            List[index] = re.sub('[A-Z]+', '', item)
    Entrez.api_key = API
    Entrez.email = email
    search_results = Entrez.read(Entrez.epost(Database, id=",".join(List)))
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    History = (webenv, query_key)

    return History

#----------------------------------------

def Find_Taxonomy_FromID(ID_List, Terms, Database, Output_File, History, email, API=None):
    from Bio import Entrez
    Batch_Size = 500
    Count = len(ID_List)
    Entrez.api_key = API
    Entrez.email = email
    Taxonomy_Dictionary = {}
    #Download in batches
    for start in range(0, Count, Batch_Size):
        end = min(Count, start+Batch_Size)
        print("Downloading record {} to {}".format(start+1, end))
        fetch_handle = Entrez.esummary(db=Database,
                                    rettype="native", retmode="text",
                                    retstart=start, retmax=Batch_Size,
                                    webenv=History[0], query_key=History[1])
        Data_Dicts = Entrez.read(fetch_handle)
        for Dict in Data_Dicts:
            for term in Terms:
            TaxID = Dict['TaxId']
            Title = Dict['Title']
            AccessionVersion = Dict['AccessionVersion']
            Taxonomy_Dictionary[AccessionVersion] = [str(TaxID), Title]
        fetch_handle.close()

    return Taxonomy_Dictionary

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script takes an NCBI ID returns some information from the database you specify.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -l [ID1, ID2, ID3] OR -f [File with IDs] -o [Output File]] -d [Database, nuccore by default]\n'''
            '''-t [Term1 Term2 Term3] -e [email] -a [API]\n'''
            '''Global mandatory parameters: -l [ID1, ID2, ID3] OR -f [File with IDs] -o [Output FastA]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-l', '--list', dest='ID_List', action='store', nargs='+', required=False, help='ID list, Separated by spaces.')
    parser.add_argument('-f', '--fileID', dest='ID_File', action='store', required=False, help='File with IDs to search, one per line')
    parser.add_argument('-d', '--db', dest='Database', action='store', required=False, default="nuccore", help='Database to look from, by default "nuccore"')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file')
    parser.add_argument('-e', '--email', dest='Email', action='store', required=True, help='Email to tell NCBI who you are.')
    parser.add_argument('-a', '--api', dest='API', action='store', required=False, help='API key for large queries')
    parser.add_argument('-t', '--terms', dest='Search_Terms', action='store', nargs='+', required=True, help='Terms to search for in the database')
    args = parser.parse_args()

    Database = args.Database.lower()
    ID_List = args.ID_List
    ID_File = args.ID_File
    Output_File = args.Output_File
    Email = args.Email.lower()
    API = args.API.lower()
    Search_Terms = args.Search_Terms

    #Parse input IDs
    Input_ID_List = []
    if ID_List == None and ID_File == None:
        sys.exit('No ID provided, either provide a list with -l or a file with IDs (one per line) with -f')
    elif ID_File != None:
        with open(ID_File) as Input_File:
            for line in Input_File:
                line = line.strip().split()[0]
                Input_ID_List.append(line)
    elif ID_List != None:
        Input_ID_List = ID_List.split(",")

    History = Post_IDs(Input_ID_List, Database, Email, API)
    Fasta_from_ID(Input_ID_List, Search_Terms, Database, Output_File, History, Email, API)

if __name__ == "__main__":
    main()