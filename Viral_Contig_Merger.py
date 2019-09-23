#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 Sept 23 2019

# Description: This script parses an MCL clustering output and returns
# the clustered and unclustered contigs
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def MCL_to_List(MCL_Output, Output_Prefix):
    """ Parses a MCL output and returns a dictionary with contigs that belong
        to each cluster and writes a FastA file with unclustered contigs. """
    import fileinput
    
    Output_dict = {}
    Unclustered_Files = []
    with open(MCL_Output) as MCL:
        Cluster = 1
        for line in MCL:
            line = line.strip().split()
            if len(line) > 1:
                if "Cluster_{}".format(Cluster) not in Output_dict:
                    Output_dict["Cluster_{}".format(Cluster)] = []
                    for genome in line:
                        Output_dict["Cluster_{}".format(Cluster)].append(genome)
                    Cluster += 1
            else:
                for genome in line:
                    Unclustered_Files.append(genome)
    print(len(Unclustered_Files))
    with open("{}_unclustered.fasta".format(Output_Prefix), 'w') as Unclust:
        for element in Unclustered_Files:
            with open(element) as input_file:
                Unclust.write(input_file.read())

    return Output_dict
    # ------------------------

def Cluster_to_Alignment(Cluster_Dict, Output_Prefix, Extension=".fa"):
    import os
    #from pymummer import coords_file, alignment, nucmer
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    #from Bio.Seq import Seq
    import pathlib

    with open("{}_Cluster_Rep.fasta".format(Output_Prefix), 'w') as Output:
        for cluster, genomes in Cluster_Dict.items():
            Sizes = []
            for genome in genomes:
                path = pathlib.Path(genome)
                with open(path, 'r') as Fasta:
                    for title, seq in SimpleFastaParser(Fasta):
                        Sizes.append((title, len(seq)))
            Sizes.sort(key=lambda tup: tup[1], reverse=True)
            Seed_genome = None
            for index, value in enumerate(Sizes):
                if index == 0:
                    genome = value[0] + Extension
                    Seed_genome = pathlib.Path(genome)
                    with open(Seed_genome) as Fasta_Genome:
                        for title, seq in SimpleFastaParser(Fasta_Genome):
                            Output.write(">{}\n{}\n".format(title,seq))

# TODO Try to dereplicate genomes...
                # else:
                #     Query_name = value[0] + ".fa"
                #     Query_genome = pathlib.Path(Contig_Folder) / Query_name
                #     runner = nucmer.Runner(Seed_genome, Query_genome, "Temp_Alignment_{}.coord".format(number), min_id = 95)
                #     runner.run()
                #     number+=1
                #     file_reader = coords_file.reader("Temp_Alignment_{}.coord".format(number))
                #     alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hit


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse
    from sys import argv

    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a MCL output into clutered and unclustered contigs\n'''
            '''Usage: ''' + argv[0] + ''' -m [MCL Output] -d [Contig Directory] -p [Output Prefix] -e [Extension of files]\n'''
            '''Global mandatory parameters: -m [MCL Output] -d [Contig Directory] -p [Output Prefix]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    parser.add_argument('-m', '--mclFile', dest='MCL_File', action='store', required=True, help='Input MCL file.')
    #parser.add_argument('-d', '--directory', dest='Contig_Dir', action='store', required=True, help='Directory where contigs are located (one per file).')
    parser.add_argument('-p', '--prefix', dest='Output_Prefix', action='store', required=True, help='Prefix of output files.')
    parser.add_argument('-e', '--extension', dest='Extension', action='store', required=True, help='Extension of contig files')
    args = parser.parse_args()

    MCL_File = args.MCL_File
    Contig_Dir = args.Contig_Dir
    Output_Prefix = args.Output_Prefix
    Extension = args.Extension

    Clusters = MCL_to_List(MCL_File, Output_Prefix)
    Cluster_to_Alignment(Clusters, Output_Prefix, Extension)

if __name__ == "__main__":
    main()