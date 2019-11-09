#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 August 4 2019

# Description: This script takes an NCBI ID of a genome and downloads the FastA file.
If you have more than 500 sequences to be downloaded, it is better to split them.
########################################################################
"""
################################################################################
"""---1.0 Import Modules---"""
from Entrez_TaxonomyLineage_fromID import Post_IDs

################################################################################
"""---2.0 Define Functions---"""

def Get_FastA_FromID(TaxID_List, History, Output, email, API = None):
    from Bio import SeqIO
    from Bio import Entrez

    Batch_Size = 500
    out_handle = open(Output, "w")
    Count = len(TaxID_List)
    Entrez.api_key = API
    Entrez.email = email
    for start in range(0, Count, Batch_Size):
        end = min(Count, start+Batch_Size)
        print("Going to download fasta record {} to {}".format(start+1, end))
        try:
            fetch_handle = Entrez.efetch(db="genome", rettype="fasta", idtype="acc", #retmode="text",
                        retstart=start, retmax=Batch_Size,
                        webenv=History[0], query_key=History[1])
            data = fetch_handle.read()
            fetch_handle.close()
            out_handle.write(data)
        except:
            print("There was an error with some of your IDs")
    out_handle.close()

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    import pandas as pd
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script takes an NCBI ID of a genome and returns the fastA file.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -l [List of IDs] or -f [File w IDs] -e [Email to tell NCBI who you are. IMPORTANT] -o [Output_Table (For list)] -a [API Key obtained from NCBI]\n'''
            '''Global mandatory parameters: -e [Email]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-l', '--list', dest='ID_List', action='store', nargs='+', required=False, help='Comma-separated list IDs, e.g., ID1,ID2,ID3...')
    parser.add_argument('-f', '--fileID', dest='ID_File', action='store', required=False, help='File with IDs to search, one per line')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='File to store all fasta sequences returned')
    parser.add_argument('-e', '--email', dest='Email', action='store', required=True, help='Email to tell NCBI who you are.')
    parser.add_argument('-a', '--api', dest='API', action='store', required=False, help='API key for large queries')
    args = parser.parse_args()

    ID_List = args.ID_List
    ID_File = args.ID_File
    Output_File = args.Output_File
    Email = args.Email.lower()
    API = args.API.lower()

    #Parse input IDs
    print(ID_List)
    if ID_File == None and ID_List == None:
        print("No IDs provided. Please give a list of IDs with -l or a file containing the IDs one per line with -f.")
    elif ID_File != None:
        Input_ID_List = []
        with open(ID_File) as Input_File:
            for line in Input_File:
                line = line.strip().split()[0]
                Input_ID_List.append(line)
    elif ID_List != None:
        Input_ID_List = ID_List

    Entry_History = Post_IDs(Input_ID_List, "assembly", Email, API)

    Get_FastA_FromID(Input_ID_List, Entry_History, Output_File, Email, API)


if __name__ == "__main__":
    main()
