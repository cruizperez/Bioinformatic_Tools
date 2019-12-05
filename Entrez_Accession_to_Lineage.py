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
def Get_UID_from_Accession(ID_Table, email, API=None):
    from time import sleep
    from Bio import Entrez
    if API != None:
        Entrez.api_key = API
    Entrez.email = email
    ID_List = list(ID_Table['Accession'])
    for Accession in ID_List:
        print(Accession)
        fetch_handle = Entrez.esearch(db='assembly', term=Accession)
        Data_Dicts = Entrez.read(fetch_handle)
        ID_Table.loc[ID_Table['Accession'] == Accession, 'UID'] = Data_Dicts["IdList"][0]
        fetch_handle.close()
        sleep(0.1)
    return ID_Table
# --------------------------------------

def Get_TaxID_from_UID(ID_Table, email, API=None):
    from time import sleep
    from Bio import Entrez
    if API != None:
        Entrez.api_key = API
    Entrez.email = email
    ID_List = list(ID_Table['UID'])
    for Accession in ID_List:
        print(Accession)
        fetch_handle = Entrez.elink(dbfrom='assembly', db="taxonomy",
                                    id=Accession)
        Data_Dicts = Entrez.read(fetch_handle)
        TaxID = Data_Dicts[0]['LinkSetDb'][0]['Link'][0]["Id"]
        ID_Table.loc[ID_Table['UID'] == Accession, 'TaxID'] = TaxID
        fetch_handle.close()
        sleep(0.1)
    return ID_Table
#---------------------------------------

def Get_Lineage_from_TaxID(ID_Table, email, API=None):
    from time import sleep
    from Bio import Entrez
    if API != None:
        Entrez.api_key = API
    Entrez.email = email
    ID_List = list(ID_Table['TaxID'])
    Valid_Ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for Accession in ID_List:
        print(Accession)
        fetch_handle = Entrez.efetch(db='taxonomy', id=Accession)
        Data_Dicts = Entrez.read(fetch_handle)
        ID_Table.loc[ID_Table['TaxID'] == Accession, ['species']] = Data_Dicts[0]['ScientificName']
        for Rank in Data_Dicts[0]['LineageEx']:
            if Rank['Rank'] in Valid_Ranks:
                ID_Table.loc[ID_Table['TaxID'] == Accession, Rank['Rank']] = Rank['ScientificName']
            else:
                pass
        fetch_handle.close()
        sleep(0.1)
    return ID_Table

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
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file')
    parser.add_argument('-e', '--email', dest='Email', action='store', required=True, help='Email to tell NCBI who you are.')
    parser.add_argument('-a', '--api', dest='API', action='store', required=False, help='API key for large queries')
    args = parser.parse_args()

    ID_List = args.ID_List
    ID_File = args.ID_File
    Output_File = args.Output_File
    Email = args.Email.lower()
    API = args.API
    if API != None:
        API = API.lower()

    # Create Empty DataFrame
    Accession_Lineage = pd.DataFrame(columns=['Accession', 'UID', 'TaxID', 'superkingdom', 'phylum', 
    'class', 'order', 'family', 'genus', 'species'])
    # Parse input IDs
    Input_ID_List = []
    if ID_List == None and ID_File == None:
        sys.exit('No ID provided, either provide a list with -l or a file with IDs (one per line) with -f')
    elif ID_File != None:
        with open(ID_File) as Input_File:
            for line in Input_File:
                line = line.strip().split()[0]
                Input_ID_List.append(line)
    elif ID_List != None:
        Input_ID_List = ID_List
    # Add IDs to DataFrame
    Accession_Lineage['Accession'] = Input_ID_List

    # Run functions
    print("Linking accessions to UIDs")
    Accession_Lineage = Get_UID_from_Accession(Accession_Lineage, Email, API)
    print("Linking UIDs to TaxIDs")
    Accession_Lineage = Get_TaxID_from_UID(Accession_Lineage, Email, API)
    print("Linking TaxIDs to Lineages")
    Accession_Lineage = Get_Lineage_from_TaxID(Accession_Lineage, Email, API)
    Accession_Lineage.to_csv(Output_File, sep="\t", header=True, index=False)

if __name__ == "__main__":
    main()