#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 August 1 2019

# Description: This script builds a sqlite database from a tab separated table.
By default it assumes the first line of the input table has headers, if not
you must provide a list of headers as -h header1,header2,header3...
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Create_Table(Database, Input_Table, Table_Name, Column_Names = None, Index_Name):
    import sqlite3
    import csv

    Header = False
    columns = None
    conn = sqlite3.connect(Database)
    conn.text_factory = str  # allows utf-8 data to be stored
    cur = conn.cursor()
    # If no columns provided extract them from the first line.
    if Column_Names == None:
        Header = True
        f = open(Input_Table)
        columns = f.readline()
        columns = columns.split("\t")
        f.close()
    else:
        columns = Column_Names
    print("Building table with columns", columns)
    # Remove previous database if exists.
    sql = "DROP TABLE IF EXISTS %s" % Table_Name
    cur.execute(sql)
    # Create new table with given column names
    sql = "CREATE TABLE %s (%s)" % (Table_Name,
                           ", ".join(["%s text" % column for column in columns]))
    cur.execute(sql)
    # Create index from the first column
    index = "%s" % (Index_Name)
    sql = "CREATE INDEX %s on %s (%s)" % ( index, Table_Name, columns[0] )
    cur.execute(sql)

    # Fill table with fields from inp
    with open(Input_Table) as DB_Table:
        reader = csv.reader(DB_Table, delimiter='\t')
        for row in reader:
            if Header == True:
                Header = False
            else:
                insertsql = "INSERT INTO %s VALUES (%s)" % (Table_Name,
                                               ", ".join([ "?" for column in row ]))
                cur.execute(insertsql, row)
        conn.commit()

    cur.close()
    conn.close()

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script builds a sqlite database from a tab separated table.\n'''
            '''By default it assumes the first line of the input table has headers, if not\n'''
            '''you must provide a list of headers as -h header1,header2,header3...\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input Table] -d [Database Name] -t [Table Name] -n [Index Name] --header [Header List]\n'''
            '''Global mandatory parameters: -i [Input Table] -d [Database Name] -t [Table Name] -n [Index Name]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Input_Table', action='store', required=True, help='Input tab-delimited table to parse, by default assumes headers are present')
    parser.add_argument('-d', '--database', dest='Database', action='store', required=True, help='Database name, can be exisiting or new.')
    parser.add_argument('-t', '--table', dest='Table_name', action='store', required=True, help='Table within db to create, always overwrites if tbale with the same name is present.')
    parser.add_argument('-n', '--index_name', dest='Index_Name', action='store', required=True, help='Name to index the table, must be unique within the database')
    parser.add_argument('--header', dest='Headers', action='store', nargs='+', required=False, help='Header list if not present in input table, enter as --header header1,header2,...')
    args = parser.parse_args()

    Input_Table = args.Input_Table
    Database = args.Database
    Table_name = args.Table_name
    Index_Name = args.Index_Name
    Headers = args.Headers

    Create_Table(Database, Input_Table, Table_name, Headers, Index_Name)


if __name__ == "__main__":
    main()
