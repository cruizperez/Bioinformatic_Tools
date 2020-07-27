#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.9
# Date:		   27 July 2020

# Description: This script parses the output from Graphbin and recovers
# the list of contigs per bin in separated files.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
import sys, argparse, os

################################################################################

"""---2.0 Define Functions---"""
def get_contig_list(input_file):
    contig_groups = {}
    with open(input_file, 'r') as bin_result:
        for line in bin_result:
            line = line.strip().split(",")
            if line[1] not in contig_groups:
                contig_groups[line[1]] = [line[0]]
            else:
                contig_groups[line[1]].append(line[0])
    
    for group, contigs in contig_groups.items():
        outfile = "MAG_" + str(group)
        with open(outfile, 'w') as output:
            for contig in contigs:
                output.write("{}\n".format(contig))

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses the output from Graphbin and recovers\n''' 
                        '''the list of contigs per bin in separated files.\n'''
                        '''Global mandatory parameters: -i [Input File]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_file', action='store', required=True, help='Graphbin output file.')
    args = parser.parse_args()

    input_file = args.input_file
    
    get_contig_list(input_file)

if __name__ == "__main__":
    main()
