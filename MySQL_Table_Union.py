
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


"""---1.0 Define Functions---"""

def Table_Inserter(Database_Name, Target_Table, List_Tables):
    # Connect to the database
    Database_Annotation = mysql.connector.connect(
                            host = "localhost",
                            user = "root",
                            passwd = "Iwillbethebest2018!",
                            database = Database_Name
                            )

    my_cursor = Database_Annotation.cursor(buffered=True)

    for Table in List_Tables:
        command = 'INSERT IGNORE INTO {}.{} (SELECT * FROM {}.{})'.format(Database_Name, Target_Table, Database_Name, Table)
        print(command)
        my_cursor.execute(command)
        Database_Annotation.commit()
        #try:
        #    my_cursor.execute(command)
        #    Database_Annotation.commit()
        #    print("Success")
        #except:
        #    print("Something went wrong")

################################################################################
"""---2.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Inserts table entries into a given table'''
                                    'Global mandatory parameters: [Matrix]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-d", "--Database_Name", dest='Database_Name', required=True, help="Name of the MySQL Database")
    parser.add_argument("-t", "--Target_Table", dest='Target_Table', required=True, help="Name of the target table inside the MySQL Database")
    parser.add_argument('-l', '--List_Tables', dest='List_Tables', required=True, nargs='+', help='List of names of the tables to insert')
    args = parser.parse_args()

    Database_Name = args.Database_Name
    Target_Table = args.Target_Table
    List_Tables = args.List_Tables

    Table_Inserter(Database_Name, Target_Table, List_Tables)

if __name__ == "__main__":
    main()
