#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz-Perez
# Email:        cruizperez3@gatech.edu
# Github:       https://github.com/cruizperez
# Institution:  Georgia Institute of Technology
# Version:      0.1
# Date:         15 February 2020

# Description: This script filters a Blast tabular output based
# on the id percentage of the read and the fraction of the match aligned (if able).
# The rapid option (--rapid) assumes you have shuffled and re-sorted the input file
# and only retains the first ocurrence per read (the highest scoring).
# The regular option compares each ocurrence of a hit and determines the best
# hit without you having to do additional work.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""

import sys, argparse, os
from random import choice
from sys import argv

################################################################################

"""---2.0 Define Functions---"""
def hit_confidence(hit, id_perc, bitscore, evalue, aln_percent, shorter, query, subject):
    if aln_percent == None:
        if float(hit[2]) >= float(id_perc) and float(hit[11]) >= float(bitscore) and float(hit[10]) <= float(evalue):
            high_quality_match = True
        else:
            high_quality_match = False
    else:
        if shorter == True:
            if (float(hit[3])*100/min(float(hit[13]),float(hit[12])) >= float(aln_percent)) \
                and float(hit[11]) >= float(bitscore) and float(hit[2]) >= float(id_perc) and float(hit[10]) <= float(evalue):
                high_quality_match = True
            else:
                high_quality_match = False
        elif query == True:
            if (float(hit[3])*100/float(hit[12])) >= float(aln_percent) and float(hit[11]) \
                >= float(bitscore) and float(hit[2]) >= float(id_perc) and float(hit[10]) <= float(evalue):
                high_quality_match = True
            else:
                high_quality_match = False
        elif subject == True:
            if (float(hit[3])*100/float(hit[13])) >= float(aln_percent) and float(hit[11]) \
                >= float(bitscore) and float(hit[2]) >= float(id_perc) and float(hit[10]) <= float(evalue):
                high_quality_match = True
            else:
                high_quality_match = False
        else:
            if (float(hit[3])*100/max(float(hit[13]),float(hit[12])) >= float(aln_percent)) \
                and float(hit[11]) >= float(bitscore) and float(hit[2]) >= float(id_perc) and float(hit[10]) <= float(evalue):
                high_quality_match = True
            else:
                high_quality_match = False
    return(high_quality_match)

### ----------------------------- Blast Parser ----------------------------------------
def blast_filter_slow(input_tab, outfile, id_perc, bitscore, evalue, aln_percent, shorter, query, subject):
    blast_hits = {}
    print("Reading " + input_tab + " Blast Output")
    # Check if Blast output has qlen and slen in addition to std output.
    with open(input_tab) as blast_input:
        if len(blast_input.readline().strip().split("\t")) == 14:
            print("I am assuming your Blast output has qlen (col 13) and slen (col 14) besides the standard output columns.")
            print("If this is not the case, i.e. columns 13 and 14 represent other values, please remove them.")
        else:
            print("I am assuming your Blast output has the standard Blast tabular output")
    # Check if any parameter was given for match filtering.
    if all(variable is None for variable in [id_perc, bitscore, evalue, aln_percent]):
        print("Retrieving best match per sequence using default parameters.")
    else:
        print("Retrieving best match per seqeunce using parameters you provided. If only some provided, the others will take their default values.")
    # Setting default values if not provided
    if id_perc is None:
        id_perc = 30
    if bitscore is None:
        bitscore = 50
    if evalue is None:
        evalue = 10
    # Retrieve best matches
    with open(input_tab) as blast_input:
        for line in blast_input:
            line = line.strip()
            hit = line.split("\t")
            good_hit = hit_confidence(hit, id_perc, bitscore, evalue, aln_percent, shorter, query, subject)
            if good_hit == True:
                if hit[0] not in blast_hits:
                    blast_hits[hit[0]] = [float(hit[11]), [line]]
                else:
                    if float(hit[11]) < blast_hits[hit[0]][0]:
                        continue
                    elif float(hit[11]) > blast_hits[hit[0]][0]:
                        blast_hits[hit[0]] = [float(hit[11]), [line]]
                    else:
                        blast_hits[hit[0]][1].append(line)
            else:
                continue
    with open(outfile, 'w') as output:
        for hit_values in blast_hits.values():
            output.write("{}\n".format(choice(hit_values[1])))
    print("Done! Check your output {}".format(outfile))

def blast_filter_fast(input_tab, outfile, id_perc, bitscore, evalue, aln_percent, shorter, query, subject):
    blast_hits = []
    print("Reading " + input_tab + " Blast Output")
    # Check if Blast output has qlen and slen in addition to std output.
    with open(input_tab) as blast_input:
        if len(blast_input.readline().strip().split("\t")) == 14:
            print("I am assuming your Blast output has qlen (col 13) and slen (col 14) besides the standard output columns.")
            print("If this is not the case, i.e. columns 13 and 14 represent other values, please remove them.")
        else:
            print("I am assuming your Blast output has the standard Blast tabular output")
    # Check if any parameter was given for match filtering.
    if all(variable is None for variable in [id_perc, bitscore, evalue, aln_percent]):
        print("Retrieving best match per sequence using default parameters.")
    else:
        print("Retrieving best match per seqeunce using parameters you provided (if only some, the others will take their default values.")
    # Setting default values if not provided
    if id_perc is None:
        id_perc = 30
    if bitscore is None:
        bitscore = 50
    if evalue is None:
        evalue = 10
    # Retrieve best matches
    with open(input_tab, 'r') as tabular, open(outfile, 'w') as output:
        for line in tabular:
            line = line.strip()
            hit = line.split("\t")
            good_hit = hit_confidence(hit, id_perc, bitscore, evalue, aln_percent, shorter, query, subject)
            if good_hit == True:
                if hit[0] not in blast_hits:
                    output.write("{}\n".format(line))
                    blast_hits.append(hit[0])
                else:
                    continue
    print("Done! Check your output {}".format(outfile))


### ------------------------------- Main function ------------------------------
def main():
        parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''Description: This script filters a Blast tabular output based\n'''
                        '''on the id percentage of the read and the fraction of the match aligned (if able).\n'''
                        '''The rapid option (--rapid) assumes you have shuffled and re-sorted the input file\n'''
                        '''and only retains the first ocurrence per read (the highest scoring).\n'''
                        '''(usually the highest scoring).\n'''
                        '''The regular option compares each ocurrence of a hit and determines the best\n'''
                        '''hit without you having to do additional work.\n'''
                        '''Basic usage: ''' + argv[0] + ''' -i [Input File] -o [Output File]\n'''
                        '''Global mandatory parameters: -i [Input File] -o [Output File]\n'''
                        '''Optional Database Parameters: See ''' + argv[0] + ' -h')
        input_options = parser.add_argument_group('Mandatory i/o options')
        input_options.add_argument('-i', '--input_tab', dest='input_tab', action='store', required=True,
                            help='Input blast result in tabular format')
        input_options.add_argument('-o', '--outfile', dest='outfile', action='store', required=True,
                            help='Output filtered blast result in tabular format')
        thresholds = parser.add_argument_group('Thresholds used for filtering')
        thresholds.add_argument('--id_perc', dest='id_perc', action='store', type=float,
                            help='Minimum percentage identity for a match to be included. By default 30')
        thresholds.add_argument('--bistcore', dest='bitscore', action='store', type=float,
                            help='Minimum bitscore for a match to be included. By default 50')
        thresholds.add_argument('--evalue', dest='evalue', action='store', type=float,
                            help='Maximum Evalue for a match to be included. By default 10')
        thresholds.add_argument('--aln_percent', dest='aln_percent', action='store', type=float,
                            help='Minimum alignment the match must cover to be included.\n' 
                                'Only calculated if you have qlen and slen in your output.')
        flags = parser.add_argument_group('Additional flags. Only calculated if you have qlen and slen in your output.')
        flags.add_argument('--longer', action='store_true',
                            help='Calculates the alignment percentage on the longer sequence, by default calculated on the shorter')
        flags.add_argument('--query', action='store_true',
                            help='Calculates the alignment percentage on the query sequence, by default false, i.e. calculated on the shorter')
        flags.add_argument('--subject', action='store_true',
                            help='Calculates the alignment percentage on the subject sequence, by default false, i.e. calculated on the shorter')
        mode = parser.add_argument_group('Filtering mode. Activate fast filtering by passing "--rapid"')
        mode.add_argument('--rapid', action='store_true',
                            help='Performs rapid filtering, see help for requirements.')
        args = parser.parse_args()
#! CORRECT LONGER VS SHORTER SEQUENCE
        input_tab = args.input_tab
        outfile = args.outfile
        id_perc = args.id_perc
        bitscore = args.bitscore
        evalue = args.evalue
        aln_percent = args.aln_percent
        longer = args.longer
        query = args.query
        subject = args.subject
        rapid = args.rapid

        if rapid == True:
            blast_filter_fast(input_tab, outfile, id_perc, bitscore, evalue, aln_percent, longer, query, subject)
        else:
            blast_filter_slow(input_tab, outfile, id_perc, bitscore, evalue, aln_percent, longer, query, subject)

if __name__ == "__main__":
    main()
