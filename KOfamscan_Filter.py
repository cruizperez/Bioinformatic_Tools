#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   14 February 2020

# Description: This script parses a KOfamscan result and returns only
# significant matches that pass the threshold.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
import argparse, sys
from random import randint

################################################################################

"""---2.0 Define Functions---"""
def kofamscan_filter(kofamscan_input, outfile):

    header_list = []
    kofamscan_results = {}
    with open(kofamscan_input, 'r') as input:
        for line in input:
            if line.startswith("#"):
                if 'gene name' in line:
                    headers = []
                    headers.append('')
                    head = line.strip().split()
                    name = ' '.join(head[0:3])
                    headers.append(name)
                    headers += head[3:]
                    header_list.append("\t".join(headers))
                else:
                    header = ['']
                    header += line.strip().split()
                    header_list.append("\t".join(header))
            elif line.startswith('*'):
                result = line.strip().split()
                gene_name = result[1]
                annot_desc = ' '.join(result[6:])
                if gene_name not in kofamscan_results:
                    kofamscan_results[gene_name] = result[0:6]
                    kofamscan_results[gene_name].append(annot_desc)
                else:
                    if result[4] > kofamscan_results[gene_name][4]:
                        kofamscan_results[gene_name] = result[0:6]
                        kofamscan_results[gene_name].append(annot_desc)
                    elif result[4] == kofamscan_results[gene_name][4]:
                        if randint(0,1) > 0:
                            kofamscan_results[gene_name] = result[0:6]
                            kofamscan_results[gene_name].append(annot_desc)
                        else:
                            continue
                    else:
                        continue
            else:
                continue

    with open(outfile, 'w') as output:
        for element in header_list:
            output.write("{}\n".format(element))
        for hit in kofamscan_results.values():
            output.write("{}\n".format("\t".join(hit)))

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a hmmsearch (HMMER) tabular ouput and\n'''
                        '''returns the best hit per gene based on the best domain OR sequence bitscore.\n'''
                        '''Global mandatory parameters: [HMMSearch Tab File] [Output File]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input", dest='kofamscan', action='store', 
                        required=True, help="Input kofamscan result")
    parser.add_argument('-o', '--output', dest='output_file', action='store', 
                        required=True, help='File to store filtered kofamscan results')
    args = parser.parse_args()

    kofamscan = args.kofamscan
    output_file = args.output_file

    # Filter kofamscan file
    kofamscan_filter(kofamscan, output_file)

if __name__ == "__main__":
    main()
