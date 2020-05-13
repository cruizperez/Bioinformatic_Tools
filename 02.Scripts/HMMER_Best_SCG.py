#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   07 February 2020

# Description: This script parses a hmmsearch (HMMER) tabular ouput and
# returns the best hit per gene per SGC.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
from random import randint

################################################################################

"""---2.0 Define Functions---"""
def hmmsearch_best_hit(hmmsearch_file, sequence, threshold, outfile):
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
                if sequence == False:
                    score = float(result[8])
                else:
                    score = float(result[5])
                if score < threshold:
                    continue
                elif result[3] in scores:
                    if score > scores[result[3]][0]:
                        scores[result[3]] = [score, line]
                    elif score == scores[result[3]][0]:
                        if randint(0,1) > 0:
                            scores[result[3]] = [score, line]
                        else:
                            continue
                    else:
                        continue
                else:
                    scores[result[3]] = [score, line]
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
            description='''This script parses a hmmsearch (HMMER) tabular ouput against a SCG dataset and\n'''
                        '''returns the best hit per SCG based on the best domain OR sequence bitscore.\n'''
                        '''Global mandatory parameters: -i [HMMSearch Tab File] -o [Output File]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--input_hmmsearch", dest='hmmsearch', action='store', 
                        required=True, help="Input hmmsearch tabular output")
    parser.add_argument('-o', '--output_file', dest='output_file', action='store', 
                        required=True, help='File to store filtered hmmsearch results')
    parser.add_argument('--sequence', dest='sequence', action='store_true', 
                        required=False, help='Filter hits by whole sequence score. Default use the domain score')
    parser.add_argument('--threshold', dest='threshold', action='store', type=int, default=0, 
                        required=False, help='Score threshold to filter. By default 0 (include all hits)')
    args = parser.parse_args()

    hmmsearch = args.hmmsearch
    output_file = args.output_file
    sequence = args.sequence
    threshold = args.threshold


    # Filter hmmsearch file
    hmmsearch_best_hit(hmmsearch, sequence, threshold, output_file)

if __name__ == "__main__":
    main()