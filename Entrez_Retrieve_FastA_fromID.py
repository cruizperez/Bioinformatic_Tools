#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 August 20 2019

# Description: This script downloads FastA records from NCBI Nucleotide DB
# from a list of accession numbers.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""
from Bio import Entrez
from Entrez_TaxonomyLineage_fromID import Post_IDs

################################################################################
"""---2.0 Define Functions---"""
def Fasta_from_ID(List, Database, Output_Fasta, History, email, API= None):
    if API != None:
        Entrez.api_key = API
    Entrez.email = email
    Batch_Size = 20
    Count = len(List)
    Output_Handle = open(Output_Fasta , 'w')
    for start in range(0, Count, Batch_Size):
        end = min(Count, start+Batch_Size)
        print("Downloading record {} to {}".format(start+1, end))
        fetch_handle = Entrez.efetch(db=Database, rettype="fasta", retmode="text",
                                     retstart=start, retmax=Batch_Size,
                                     webenv=History[0], query_key=History[1])
        Data = fetch_handle.read()
        fetch_handle.close()
        Output_Handle.write(Data)
    Output_Handle.close()

#--------------------------------------

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script takes an NCBI ID of a genome and returns the FastA file.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -l [ID1, ID2, ID3] -o [Output FastA] -d [Database, nuccore by default]\n'''
            '''Global mandatory parameters: -l [ID1, ID2, ID3] -o [Output FastA]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-l', '--list', dest='ID_List', action='store', nargs='+', required=True, help='ID list, Separated by spaces.')
    parser.add_argument('-d', '--db', dest='Database', action='store', required=False, default="nuccore", help='Database to look from, by default "nuccore"')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file')
    parser.add_argument('-e', '--email', dest='Email', action='store', required=True, help='Email to tell NCBI who you are.')
    parser.add_argument('-a', '--api', dest='API', action='store', required=False, help='API key for large queries')
    args = parser.parse_args()

    Database = args.Database.lower()
    ID_List = args.ID_List
    Output_File = args.Output_File
    Email = args.Email.lower()
    API = args.API.lower()


    History = Post_IDs(ID_List, Database, Email, API)
    Fasta_from_ID(ID_List, Database, Output_File, History, Email, API)

if __name__ == "__main__":
    main()

