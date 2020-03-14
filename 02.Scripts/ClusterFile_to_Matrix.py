#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.1
# Date:		   07 February 2020

# Description: This script calculates the completeness and contamination
# of genomes from a filtered HMM search results file. 
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""
import argparse, sys
import pandas as pd

################################################################################
"""---2.0 Define Functions---"""

def parse_cluster_file(cluster_file, genome_separator, outfile):
    total_clusters = []
    total_genomes = []
    cluster = None
    cluster_genes = {}
    cluster_representative = {}
    with open(cluster_file) as input_file:
        for line in input_file:
            if line.startswith(">"):
                cluster = line.strip().replace(">", "")
                total_clusters.append(cluster)
                cluster_genes[cluster] = []
            else:
                line = line.strip().replace(">", "")
                line = line.replace("...", "")
                gene = line.split()[2]
                if line.split()[0] == "0":
                    cluster_representative[cluster] = gene
                elif line.split()[-1] == "*":
                    cluster_representative[cluster] = gene
                genome = genome_separator.join(gene.split(genome_separator)[0:-1])
                total_genomes.append(genome)
                cluster_genes[cluster] .append(genome)

    total_genomes = sorted(list(set(total_genomes)))
    cluster_table = pd.DataFrame(0, index=total_clusters, columns=total_genomes)
    with open(outfile+".rep", 'w') as output:
        for cluster, gene in cluster_representative.items():
            output.write("{}\t{}\n".format(cluster,gene))
    for cluster, genomes in cluster_genes.items():
        for genome in genomes:
            cluster_table.loc[cluster, genome] += 1
    cluster_table.to_csv(outfile+".matrix", sep="\t", header=True, index=True)

################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''Parses a CD-HIT-type cluster file and returns a matrix for cluster presence/absence\n'''
                        '''Global mandatory parameters: [Cluster File] [Output File]\n'''
                        '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='cluster_file', action='store', required=True, help='Cluster file to parse')
    parser.add_argument('-o', '--output', dest='output_file', action='store', required=True, help='Output prefix for the matrix counts and representative list')
    parser.add_argument('--separator', dest='separator', action='store', required=False, default="--",
                        help='String separating genome name from contig and gene, by default "--".\nIf separator has "-" or "--" pass it as --separator="--"')
    args = parser.parse_args()

    cluster_file = args.cluster_file
    output_file = args.output_file
    separator = args.separator

    if isinstance(separator, list):
        separator = separator[0]
    else:
        separator =  separator

    parse_cluster_file(cluster_file, separator, output_file)

if __name__ == "__main__":
    main()
