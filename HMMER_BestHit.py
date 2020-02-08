#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   07 February 2020

# Description: This script parses a hmmsearch (HMMER) tabular ouput and
# returns the best hit per gene based on the best domain OR sequence bitscore. 
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
from random import randint

################################################################################

"""---2.0 Define Functions---"""
def hmmsearch_best_hit(hmmsearch_file, domain, outfile):
    """
    Filters a hmmsearch tab output returning the best hit per protein
    
    Arguments:
        hmmsearch_file {filepath} -- File with hmmsearch results
    """
    headers = []
    scores = {}
    with open(hmmsearch_file, 'r') as hmm_input:
        for line in hmm_input:
            line = line.strip()
            if line.startswith("#"):
                headers.append(line)
            else:
                result = line.strip().split()
                if domain == True:
                    score = result[8]
                else:
                    score = result[5]
                if score in scores:
                    if score > scores[result[0]][0]:
                        scores[result[0]] = [score, line]
                    elif score == scores[result[0]][0]:
                        if randint(0,1) > 0:
                            scores[result[0]] = [score, line]
                        else:
                            continue
                    else:
                        continue
                else:
                    scores[result[0]] = [score, line]
    with open(outfile, 'w') as output:
        for element in headers[0:3]:
            output.write("{}\n".format(element))
        for value in scores.values():
            output.write("{}\n".format(value[1]))
        for element in headers[3:]:
            output.write("{}\n".format(element))

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a hmmsearch (HMMER) tabular ouput and\n'''
                        '''returns the best hit per gene based on the best domain OR sequence bitscore.\n'''
                        '''Global mandatory parameters: [HMMSearch Tab File] [Output File]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input_hmmsearch", dest='hmmsearch', action='store', 
                        required=True, help="Input hmmsearch tabular output")
    parser.add_argument('-o', '--output_file', dest='output_file', action='store', 
                        required=True, help='File to store filtered hmmsearch results')
    parser.add_argument('--domain', dest='domain', action='store_false', 
                        required=False, help='Filter hits by domain score, by default True')
    args = parser.parse_args()

    hmmsearch = args.hmmsearch
    output_file = args.output_file
    domain = args.domain


    # Filter hmmsearch file
    hmmsearch_best_hit(hmmsearch, domain, output_file)

if __name__ == "__main__":
    main()
