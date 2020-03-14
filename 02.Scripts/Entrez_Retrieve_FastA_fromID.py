#!/usr/bin/env python

"""
########################################################################
# Author:      Carlos Ruiz, cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# https://github.com/cruizperez/
# Version:    2.0
# Date:      07 Nov 2019

# Description: This script downloads FastA records from a given NCBI DB
# from a list of accession numbers.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""
def Fetch_FastA(List, Database, Single_Output, email=None, API=None):
    import subprocess
    if email == None:
        print("It is recommended to provide an email address (option -e).")
        email = "user@computer.com"
    if API == None:
        print("If you do not provide an API you are limited to 3 requests per second.")
        print("Running with no API provided.")
        for ID in List:
            Entrez = subprocess.Popen(["esearch", "-email", email, "-db", Database, "-query", ID, "|", "efetch", 
                                    "-format", "fasta"], stdout=subprocess.PIPE)
            Entrez.wait()
            out, err = Entrez.communicate()
            print(out)
    else:
        if Single_Output != None:
            for ID in List:
                Command = 'esearch -email {} -api {} -db {} -query {} | efetch -format fasta > {}'.format(email,
                            API, Database, ID, Single_Output)
                print(Command)
                Efetch = subprocess.Popen(Command,shell=True)#,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                Efetch.wait()
                # out= Efetch.communicate()
                # print(out)
                # Esearch_Process = subprocess.Popen(["esearch", "-email", email, "-api", API, "-db", Database, "-query", 
                                        # ID], stdout=subprocess.PIPE)
                # Entrez.wait()
                # out= Esearch_Process.communicate()
                # Efetch_Process = subprocess.Popen(["efetch", "-format", "fasta"],  stdout=subprocess.PIPE)
                # print(out)


# def Post_IDs(List, Database, email, API= None):
#     if Database == "biosample":
#         import re
#         for index, item in enumerate(List):
#             List[index] = re.sub('[A-Z]+', '', item)
#     Entrez.api_key = API
#     Entrez.email = email
#     search_results = Entrez.read(Entrez.epost(Database, id=",".join(List)))
#     webenv = search_results["WebEnv"]
#     query_key = search_results["QueryKey"]
#     History = (webenv, query_key)

#     return History

# def Fasta_from_ID(List, Database, Output_Fasta, History, email, API= None):
#     if API != None:
#         Entrez.api_key = API
#     Entrez.email = email
#     Batch_Size = 20
#     Count = len(List)
#     Output_Handle = open(Output_Fasta , 'w')
#     for start in range(0, Count, Batch_Size):
#         end = min(Count, start+Batch_Size)
#         print("Downloading record {} to {}".format(start+1, end))
#         fetch_handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text",
#                                      retstart=start, retmax=Batch_Size,
#                                      webenv=History[0], query_key=History[1])
#         Data = fetch_handle.read()
#         fetch_handle.close()
#         Output_Handle.write(Data)
#     Output_Handle.close()

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
    parser.add_argument('-l', '--list', dest='ID_List', action='store', nargs='+', required=False, help='ID list, Separated by spaces.')
    parser.add_argument('-f', '--fileID', dest='ID_File', action='store', required=False, help='File with IDs to search, one per line')
    parser.add_argument('-d', '--db', dest='Database', action='store', required=False, default="nuccore", help='Database to look from, by default "nuccore"')
    # parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file')
    parser.add_argument('-e', '--email', dest='Email', action='store', required=True, help='Email to tell NCBI who you are.')
    parser.add_argument('-a', '--api', dest='API', action='store', required=False, help='API key for large queries')
    parser.add_argument('-s', '--single', dest='Single_Output', action='store', required=False, help='Output file if you need all downloads in a single file, by default it downloads individually')

    args = parser.parse_args()

    Database = args.Database.lower()
    ID_List = args.ID_List
    ID_File = args.ID_File
    # Output_File = args.Output_File
    Email = args.Email.lower()
    API = args.API.lower()
    Single_Output = args.Single_Output

    #Parse input IDs
    if ID_List == None and ID_File == None:
        sys.exit('No ID provided, either provide a list with -l or a file with IDs (one per line) with -f')
    elif ID_File != None:
        Input_ID_List = []
        with open(ID_File) as Input_File:
            for line in Input_File:
                line = line.strip().split()[0]
                Input_ID_List.append(line)
    elif ID_List != None:
        Input_ID_List = ID_List

    Fetch_FastA(Input_ID_List, Database, Single_Output, email=Email, API=API)
    # History = Post_IDs(Input_ID_List, Database, Email, API)
    # Fasta_from_ID(Input_ID_List, Database, Output_File, History, Email, API)

if __name__ == "__main__":
    main()

