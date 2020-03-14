#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   1.0
# Date:		   11 February 2020

# Description: This script looks at a list of files (one per line)and checks
# if each file exists in a given location. Then, it creates a second list 
# with the unexisting files.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
from pathlib import Path

################################################################################

"""---2.0 Define Functions---"""
def check_files_exist(input_list, folder, outfile):
    with open(input_list, 'r') as input, open(outfile, 'w') as output:
        folder_path = Path(folder)
        for line in input:
            file_id = line.strip().split("\t")[0]
            complete_path = folder_path / file_id
            if complete_path.exists():
                continue
            else:
                output.write("{}\n".format(file_id))

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script looks at a list of files (one per line)and checks\n'''
                        '''if each file exists in a given location. Then, it creates a second list\n'''
                        '''with the unexisting files.\n'''
                        '''Global mandatory parameters: -i [File list] -f [Folder to search] -o [Output file]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input_list", dest='input_list', action='store', 
                        required=True, help="Input list with file names")
    parser.add_argument('-f', '--folder', dest='folder', action='store', 
                        required=True, help='Folder to search for the files')
    parser.add_argument('-o', '--output', dest='outfile', action='store', 
                        required=False, help='File to save list of unexisting files')
    args = parser.parse_args()

    input_list = args.input_list
    folder = args.folder
    outfile = args.outfile

    check_files_exist(input_list, folder, outfile)

if __name__ == "__main__":
    main()
