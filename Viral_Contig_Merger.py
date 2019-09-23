#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  2.0
# Date:		 August 12 2019

# Description: This searches IDs provided in a given database.
# Its use is mainly to search for RefSeq and Uniprot IDs faster than parsing files.
########################################################################
"""
# ----0.0 Prerequisites-----




################################################################################
"""---1.0 Define Functions---"""

def MCL_to_List(MCL_Output):
    """ Parses a MCL output and returns a dictionary with contigs that belong
        to each cluster and writes a FastA file with unclustered contigs. """
    import fileinput
    
    Output_dict = {}
    Unclustered_Files = []
    with open("out.OV_Genomes-vs-OV_Genomes.ANI.mci.I20") as MCL:
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
    Unclustered_Files = ["01.Final_OV_Contgs/" + s for s in Unclustered_Files]
    print(len(Unclustered_Files))
    with open("Unclustered2.fasta", 'w') as Unclust:
        for element in Unclustered_Files:
            with open(element) as input_file:
                Unclust.write(input_file.read())

    return Output_dict
    # ------------------------

def Cluster_to_Alignment(MCL_Output, Contig_Folder, Output_Folder):
    import os
    import pymummer


Clusters = MCL_to_List("Input")
