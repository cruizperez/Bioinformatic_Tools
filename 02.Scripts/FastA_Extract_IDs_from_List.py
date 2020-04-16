#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Github:       https://github.com/cruizperez
# Institution:   Georgia Institute of Technology
# Version:      0.9
# Date:         March 17, 2020

# Description: This script reads a list of incomplete fasta ids from a list
# and extracts the complete id from the original fasta file.
########################################################################
"""
################################################################################
"""---0.0 Import Modules---"""
import sys, argparse

################################################################################
"""---1.0 Define Functions---"""

def extract_complete_ids(id_list, fasta_file, outfile):
    fasta_ids = []
    with open(fasta_file, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                line = line.strip().split()
                sequence_id = line[0].replace(">", "")
                fasta_ids.append(sequence_id)
    with open(id_list, 'r') as input_list, open(outfile, 'w') as output:
        for line in input_list:
            line = line.strip()
            for sequence_id in fasta_ids:
                if line in sequence_id:
                    output.write("{}\n".format(sequence_id))
                    break

################################################################################
"""---3.0 Main Function---"""

def main():
    # Description: 
# 
# 
# 
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a table with protein IDs and searches\n'''
            '''these IDs in a SQL database to extract protein annotations.\n'''
            '''You can also include a column with the ids of your proteins (queries)\n'''
            '''so they will be included in the output (use query_col [int]).\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input Table] -d [Database Name] -t [Table Name] --query_col [Query Column]\n'''
            '''Global mandatory parameters: -i [Input Table] -d [Database Name] -t [Table Name] --query_col [Query Column]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_file', action='store', required=True,
                        help='Input tab-delimited table to parse, by default assumes no headers are present')
    parser.add_argument('-o', '--output', dest='outfile', action='store', required=True,
                        help='Output tab-delimited table to store annotations')
    parser.add_argument('-d', '--database', dest='database', action='store', required=True,
                        help='SQL database where the annotations are stored')
    parser.add_argument('-t', '--table', dest='table', action='store', required=True,
                        help='Table within SQLite database to search in. Choose between "swissprot", "trembl" or "refseq"')
    parser.add_argument('--target_col', dest='target_col', action='store', required=False, type=int, default=1, 
                        help='Column with ids of database hits (targets) in the input file. By default 1, first column.')
    parser.add_argument('--query_col', dest='query_col', action='store', required=False, type=int,
                        help='Column with query ids in the input file, by default None, i.e. assumes only target ids are given.')
    args = parser.parse_args()

    input_file = args.input_file
    outfile = args.outfile
    database = args.database
    table = args.table
    target_col = args.target_col - 1
    query_col = args.query_col - 1

    # Read input table and extract ids.
    input_list = parse_input_table(input_file, query_col, target_col)
    search_table(database, table, input_list, outfile)

if __name__ == "__main__":
    main()
