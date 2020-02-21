#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   21 February 2020

# Description: This script takes a list of files and an absolute folder location
# and creates a Qiime2 manifest file.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
import argparse, sys
from pathlib import Path

################################################################################
"""---2.0 Define Functions---"""

def qiime2_manifest_creator(folder, file_list, field_separator, outfile, sample_col, read_col):
    read_files = []
    folder_path = Path(folder)
    for file in file_list:
        file = Path(file).name
        read_files.append(file)
    
    with open(outfile, 'w') as output:
        output.write("sample-id,absolute-filepath,direction\n")
        for file in read_files:
            file_parsed = file.split(field_separator)
            sample = file_parsed[sample_col-1]
            orientation = file_parsed[read_col-1]
            if orientation == "R1":
                orientation  = "forward"
            elif orientation == "R2":
                orientation = "reverse"
                
            output.write("{},{},{}\n".format(sample, str(folder_path / Path(file)), orientation))


################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script takes a list of files and an absolute folder location and creates a Qiime2 manifest file.\n'''
                        '''Global mandatory parameters: -f [Folder] -o [Output File] -i OR -l [Input files]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input_list', dest='input_list', action='store', required=False, help='File with list of files to include')
    parser.add_argument('-l', '--list_files', dest='file_list', action='store', required=False, nargs="+", help='List of files to include')
    parser.add_argument('-f', '--folder', dest='folder', action='store', required=True, help='Absolute path to folder where files are located')
    parser.add_argument('-o', '--output', dest='output_file', action='store', required=True, help='Output file to store manifest')
    parser.add_argument('--separator', dest='separator', action='store', required=False, default="_", help='String separating the fields in the file name. By default "_"')
    parser.add_argument('--sample_col', dest='sample_col', action='store', default=1, type=int, required=False,
                        help='''Column in the filename where the sample name is. By default 1.\n
                                For example, in the file name HV0627_S25_L001_R1_001.fastq separated by "_", the sample col is 1''')
    parser.add_argument('--read_col', dest='read_col', action='store', default=4, type=int, required=False,
                        help='''Column in the filename where the read orientation R1/R2 is. By default 4.\n
                                For example, in the file name HV0627_S25_L001_R1_001.fastq separated by "_", the orientation col is 4''')
    args = parser.parse_args()

    input_list = args.input_list
    file_list = args.file_list
    folder = args.folder
    output_file = args.output_file
    separator = args.separator
    sample_col = args.sample_col
    read_col = args.read_col

    if input_list != None and file_list != None:
        sys.exit("Please provide only a list of files e.g. $(ls *.fastq) OR a file with the file names to include")
    elif file_list != None:
        qiime2_manifest_creator(folder, file_list, separator, output_file, sample_col, read_col)
    elif input_list != None:
        file_list = []
        with open(input_list, 'r') as input:
            for line in input:
                line = line.strip()
                file_list.append(line)
        qiime2_manifest_creator(folder, file_list, separator, output_file, sample_col, read_col)
    else:
        sys.exit("No input provided, please use -i or -l to give me the input :)")

if __name__ == "__main__":
    main()