#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Institution:   Georgia Institute of Technology
# Version:	  1.2
# Date:		 09 February 2020

# Description: This script filters a MagicBlast tabular output based
# on the id percentage of the read and the length aligned.
# The rapid option (--rapid) assumes you have shuffled and re-sorted the input file
# and only retains the first ocurrence per read (usually the highest scoring).
# The regular option compares each ocurrence of a hit and determines the best
# hit without you having to do additional work.
########################################################################
"""
################################################################################
"""---0.0 Import Modules---"""
from numpy.random import randint
import argparse
from sys import argv

################################################################################
"""---1.0 Define Functions---"""

def MagicBlast_filter_slow(input_tab, outfile, aln_fraction = 80, percent_id = 1):
    MagicBlast_hits = {}
    with open(input_tab) as tabular:
        for line in tabular:
            line = line.strip()
            hit = line.split()
            if hit[0] == '#':
                continue
            elif float(hit[2]) < percent_id:
                continue
            elif (float(hit[7]) - float(hit[6]) * 100 / float(hit[15])) < aln_fraction:
                continue
            else:
                score = float(hit[12])
                if hit[0] not in MagicBlast_hits:
                    MagicBlast_hits[hit[0]] = [score, line]
                else:
                    if score < MagicBlast_hits[hit[0]][0]:
                        continue
                    elif score > MagicBlast_hits[hit[0]][0]:
                        MagicBlast_hits[hit[0]] = [score, line]
                    else:
                        if randint(2) > 0:
                            MagicBlast_hits[hit[0]] = [score, line]
                        else:
                            continue
    with open(outfile, 'w') as output:
        for hit_values in MagicBlast_hits.values():
            output.write("{}\n".format(hit_values[1]))
                    
def MagicBlast_filter_rapid(input_tab, outfile, aln_fraction = 80, percent_id = 1):
    MagicBlast_hits = []
    with open(input_tab, 'r') as tabular, open(outfile, 'w') as output:
        for line in tabular:
            line = line.strip()
            hit = line.split()
            if hit[0] == '#':
                continue
            elif float(hit[2]) < percent_id:
                continue
            elif (float(hit[7]) - float(hit[6]) * 100 / float(hit[15])) < aln_fraction:
                continue
            else:
                if hit[0] not in MagicBlast_hits:
                    output.write("{}\n".format(line))
                    MagicBlast_hits.append(hit[0])
                else:
                    continue

################################################################################
"""---2.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''Description: This script filters a MagicBlast tabular output based\n'''
                        '''on the id percentage of the read and the length aligned.\n'''
                        '''The rapid option (--rapid) assumes you have shuffled and re-sorted\n'''
                        '''the input file and only retains the first ocurrence per read\n'''
                        '''(usually the highest scoring).\n'''
                        '''The regular option compares each ocurrence of a hit and determines\n'''
                        '''the best hit without you having to do additional work.\n'''
            '''Usage: ''' + argv[0] + ''' -i [Input File] -o [Output File] -p [Identity] -f [Fraction Aligned]\n'''
            '''Global mandatory parameters: -i [Input File] -o [Output File]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    parser.add_argument('-i', '--input_tab', dest='input_tab', action='store', required=True, 
    help='Input MagicBlast tabular file.')
    parser.add_argument('-o', '--output', dest='outfile', action='store', required=True, 
    help='Output Table')
    parser.add_argument('-p', '--pidentity', dest='percent_id', action='store', required=True, type=int, 
    help='Percentage identity of matches to retain. By default no identity filter.', default=0)
    parser.add_argument('-f', '--fraction', dest='aln_fraction', action='store', required=False, type=int, 
    help='Minimum percentage of read aligned to be included. By default 80', default=80)
    parser.add_argument('--rapid', dest='rapid_filter', action='store_true', required=False, 
    help='Perform rapid filter. Only retains first high quality occurrence. Useful if pre-shuffled and sorted input.')
    args = parser.parse_args()

    input_tab = args.input_tab
    outfile = args.outfile
    percent_id = args.percent_id
    aln_fraction = args.aln_fraction
    rapid_filter = args.rapid_filter

    if rapid_filter == True:
        MagicBlast_filter_rapid(input_tab, outfile, aln_fraction, percent_id)
    else:
        MagicBlast_filter_slow(input_tab, outfile, aln_fraction, percent_id)

if __name__ == "__main__":
    main()