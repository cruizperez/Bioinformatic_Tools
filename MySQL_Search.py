#!/usr/bin/env python

"""
########################################################################
# Author:      Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:    1.0
# Date:      04 May 2019

# Description: This script creates a MySQL database from the NCBI bacterial,
# archaeal and viral RefSeq databases for annotation purposes.

########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import os, sys, argparse
import mysql.connector
from Bio import SeqIO


################################################################################
"""---1.0 Define Functions---"""
# Connect to the database
NR_Database_Annotation = mysql.connector.connect(
                        host = "localhost",
                        user = "root",
                        passwd = "Iwillbethebest2018!",
                        database = "ncbi_refseq"
                        )

# Create cursor to execute commands into the database
my_cursor = NR_Database_Annotation.cursor(buffered=True)


def Database_Search(Input_ID_Table):
    with open(Input_ID_Table) as ID_Table:
        for Line in ID_Table:
            Prot_Record = (line.split("\t")[1])
            print(Prot_Record)
            try:
                my_cursor.execute("""SELECT * FROM ncbi_refseq.merged WHERE MATCH(ID2) AGAINST('%s' IN NATURAL LANGUAGE MODE);""", (Prot_Record,))
            except:
                print("ID not found in DB")

################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''This script creates a MySQL database from the NCBI bacterial,\n
                                                    archaeal and viral RefSeq databases for annotation purposes.'''
                                    'Global mandatory parameters: [Genbank_File(s)] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--inputIDs', dest='Input_ID_Table', required=True, help='Blast Output File')
    #parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    args = parser.parse_args()

    Input_ID_Table = args.Input_ID_Table
    #Output_File = args.Output_File

    Database_Search(Input_ID_Table)

if __name__ == "__main__":
    main()
