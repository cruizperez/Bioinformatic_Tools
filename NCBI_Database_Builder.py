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
#from subprocess import Popen, PIPE
#import subprocess
#import xml.etree.ElementTree as ET
#from Bio import Entrez
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

## If database not created uncomment and create it.
#my_cursor.execute("CREATE DATABASE ncbi_refseq")
#my_cursor.execute("SHOW DATABASES")
#my_cursor.execute("""CREATE TABLE Annotation(
#                    ID VARCHAR(255) PRIMARY KEY,
#                    Product VARCHAR(255),
#                    Taxonomy TEXT,
#                    Note TEXT)""")


def Database_Filler(Input_Genbank):
    for File in Input_Genbank:
        print(File)
        print("----------------------------vim ")
        GBFile = SeqIO.parse(File, "genbank")
        for record in GBFile:
            ID = record.id
            print(ID)
            if 'comment' not in record.annotations:
                Note == ""
                #print(Note)
            else:
                Note = record.annotations['comment']
                Note = Note.replace("\n", " ")
                #print(Note)
            if 'organism' not in record.annotations:
                genus = 'NA'
                species = 'NA'
                #print(genus,species)
            else:
                if len(record.annotations['organism'].split()) > 1:
                    genus = record.annotations['organism'].split()[0]
                    species = record.annotations['organism']
                    #print(genus,species)
                else:
                    genus = record.annotations['organism']
                    species = ""
                    #print(genus,species)
            if 'taxonomy' not in record.annotations:
                #print("No taxonomy found")
                taxonomy = 'NA'
                #domain = 'NA'
                #phylum = 'NA'
                #class_t = 'NA'
                #order = 'NA'
                #family = 'NA'
            else:
                taxonomy = ", ".join(record.annotations['taxonomy']) + ", " + species
                #for level in range(0, len(record.annotations['taxonomy'])):
                #    if level == 0:
                #        domain = record.annotations['taxonomy'][0]
                #    elif level == 1:
                #        phylum = record.annotations['taxonomy'][1]
                #    elif level == 2:
                #        class_t = record.annotations['taxonomy'][2]
                #    elif level == 3:
                #        order = record.annotations['taxonomy'][3]
                #    elif level == 4:
                #        family = record.annotations['taxonomy'][4]
                #    else:
                #        continue
            #print(taxonomy)
            for feature in record.features:
                if feature.type == "Protein":
                    if 'product' not in feature.qualifiers:
                        if 'name' in feature.qualifiers:
                            Product = feature.qualifiers['name'][0]
                        else:
                            Product = 'NA'
                        #print(Product)
                    else:
                        Product = feature.qualifiers['product'][0]
                        #print(Product)

            Insert_mySQL = """INSERT INTO Annotation (ID, Product, Taxonomy, Note)
                                                    VALUES (%s, %s, %s, %s)"""
            #print(Product, taxonomy, Note)
            Prot_Record = (ID, Product, taxonomy, Note)
            try:
                my_cursor.execute(Insert_mySQL, Prot_Record)
                NR_Database_Annotation.commit()
            except:
                print("ID already in DB")

################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''This script creates a MySQL database from the NCBI bacterial,\n
                                                    archaeal and viral RefSeq databases for annotation purposes.'''
                                    'Global mandatory parameters: [Genbank_File(s)] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--gbfile', dest='GenbankFiles', nargs='+', required=True, help='Genabnk files to extract information from')
    #parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    args = parser.parse_args()

    GenbankFiles = args.GenbankFiles
    #Output_File = args.Output_File

    Database_Filler(GenbankFiles)

if __name__ == "__main__":
    main()
