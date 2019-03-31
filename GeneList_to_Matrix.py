#!/usr/bin/env python

"""
########################################################################
# Author:      Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:    1.0
# Date:      20 January 2019

# Description: This script parses one or more lists of genes present in a genome
# and returns a matrix with gene abundances per genome.

########################################################################
"""

"""--------------------0.0 Import Modules --------------------------------"""
import pandas as pd
import sys, argparse

"""--------------------1.0 Initialize Variables---------------------------"""

parser = argparse.ArgumentParser(description=
    'Global mandatory parameters: -i [input_file(s)] -o [output]\n')


parser.add_argument("-i", "--inputFiles", required=True, nargs='+', help="Input list of genes present per genome")
parser.add_argument("-o", "--outputFile", required=False, default="GenomeMatrix", help="Output genome matrix")

args = parser.parse_args()

inputFile = args.inputFiles
outputFile = args.outputFile

"""--------------------2.0 Analyze Files-----------------------------------"""

# Create function to interate line per line and create counts in a dictionary, then transform the dictionary into a dataframe
# Maybe I can save one step filling a dataframe right away...
def TableParser(inputFile):
    inputFile_geneMatrix = {}
    table = open(inputFile, "r")
    for line in table:
        line = line.strip()
        inputFile_geneMatrix[line] = inputFile_geneMatrix.get(line, 0) + 1
    GenomeMatrix[inputFile] = pd.DataFrame.from_dict(inputFile_geneMatrix, orient='index', columns=[inputFile])


GenomeMatrix = {}

# Iterate through each input file.
for file in inputFile:
    TableParser(file)

# Concatenate all dataframes from the dictionary into a single dataframe and export it.
result = pd.concat(GenomeMatrix.values(), axis=1, sort=False).fillna(0)
result.index.name = "Genes"
result.to_csv(outputFile, sep='\t')