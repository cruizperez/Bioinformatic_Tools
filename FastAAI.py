#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 Nov 27 2019

# Description: Calculates the average amino acid identity using 
# a tab-formatted Sword/Blast/Diamond results from a previously 
# calculated all-vs-all comparison.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Genes_Per_Genome(Input, Gene_Separator,Contig_Separator):
    """Counts the number of genes per genome and the length per gene
    of the protein file provided"""
    Number_Genes = {}
    Gene_Length = {}
    with open(Input) as FastAInput:
        for line in FastAInput:
            if ">" in line:
                Gene = line.split()[0].replace(">","")
                Gene_Length[Gene] = 0
                Genome = Gene.split(Contig_Separator)
                Genome = Contig_Separator.join(Genome[:-1])
                Number_Genes[Genome] = Number_Genes.get(Genome, 0) + 1
            else:
                line = line.strip()
                Gene_Length[Gene] += len(line)
    return (Number_Genes, Gene_Length)

def BH_per_Genome(Input, Output, Gene_Length, Gene_Separator, Contig_Separator):
    from numpy import random
    
    Best_Hit_per_Genome = {}
    with open(Input) as SearchInput:
        for line in SearchInput:
            line = line.strip().split()
            Query_Gene = line[0]
            if float(line[2]) >= 30 and (float(line[3]) / Gene_Length[Query_Gene]) >= 0.7: # Set ID and Fraction
                Reference_Gene = line[1]
                # Find the Reference Genome
                if Contig_Separator == None:
                    Ref_Genome = line[1].split(Gene_Separator)
                    Ref_Genome = Gene_Separator.join(Ref_Genome[:-1])
                else:
                    Ref_Genome = line[1].split(Contig_Separator)
                    Ref_Genome = Contig_Separator.join(Ref_Genome[:-1])

                # Find best hit per gene in each genome
                if (Query_Gene, Ref_Genome) not in Best_Hit_per_Genome:
                    Best_Hit_per_Genome[(Query_Gene, Ref_Genome)] = line
                else:
                    if float(line[11]) > float(Best_Hit_per_Genome[(Query_Gene, Ref_Genome)][11]):
                        Best_Hit_per_Genome[(Query_Gene, Ref_Genome)] = line
                    elif float(line[11]) == float(Best_Hit_per_Genome[(Query_Gene, Ref_Genome)][11]):
                        if Query_Gene == Reference_Gene:
                            Best_Hit_per_Genome[(Query_Gene, Ref_Genome)] = line
                        elif random.choice([0,1]) > 0:
                            Best_Hit_per_Genome[(Query_Gene, Ref_Genome)] = line
                        else:
                            continue
                    else:
                        continue
            else:
                continue
    # Save table with best hit per gene per genome
    with open(Output, 'w') as Output:
        for Hit in Best_Hit_per_Genome.values():
            Output.write("{}\n".format("\t".join(Hit)))


def TwoWay_AAI(Input, Output, Gene_Separator, Contig_Separator):
    from numpy import mean
    from numpy import std
    import pandas as pd

    Hit_Dictionary = {}
    # Store hits in a dictionary
    print("Parsing gene hits")
    with open(Input, 'r') as BestHits:
        for line in BestHits:
            line = line.strip().split()
            Query_Gene = line[0]
            Ref_Gene = line[1]
            if (Query_Gene, Ref_Gene) in Hit_Dictionary:
                Hit_Dictionary[(Query_Gene, Ref_Gene)].append(float(line[2]))
            elif (Ref_Gene, Query_Gene) in Hit_Dictionary:
                Hit_Dictionary[(Ref_Gene, Query_Gene)].append(float(line[2]))
            else:
                Hit_Dictionary[(Query_Gene, Ref_Gene)] = [float(line[2])]

    # Remove entries where no RBM was identified
    Key_List = list(Hit_Dictionary.keys())
    for i in Key_List:
        if i[0] == i[1]:
            continue
        elif len(Hit_Dictionary[i]) < 2:
            del Hit_Dictionary[i]

    # Parse reduced dictionary and calculate AAI
    print("Calculating two-way AAI")
    Key_List = list(Hit_Dictionary.keys())
    AAI = {}
    for item in Key_List:
        if item[0] == item[1] or len(Hit_Dictionary[item]) > 1:
            Query = item[0]
            Reference = item[1]
            if Contig_Separator == None:
                Query_Genome = Query.split(Gene_Separator)
                Query_Genome = Gene_Separator.join(Query_Genome[:-1])
                Ref_Genome = Reference.split(Gene_Separator)
                Ref_Genome = Gene_Separator.join(Ref_Genome[:-1])
            else:
                Query_Genome = Query.split(Contig_Separator)
                Query_Genome = Contig_Separator.join(Query_Genome[:-1])
                Ref_Genome = Reference.split(Contig_Separator)
                Ref_Genome = Contig_Separator.join(Ref_Genome[:-1])
            if (Query_Genome, Ref_Genome) in AAI:
                AAI[(Query_Genome, Ref_Genome)].append(mean(Hit_Dictionary[(Query, Reference)]))
            elif (Ref_Genome, Query_Genome) in AAI:
                AAI[(Ref_Genome, Query_Genome)].append(mean(Hit_Dictionary[(Query, Reference)]))
            else:
                AAI[(Query_Genome, Ref_Genome)] = [mean(Hit_Dictionary[(Query, Reference)])]

    print("Creating AAI Table")
    Hit_List = []
    for genomes, hits in AAI.items():
        Hit_List.append([genomes[0], genomes[1], round(mean(hits),2), round(std(hits),2), len(hits), min(Number_Genes[genomes[0]], Number_Genes[genomes[1]])])
    AAI_Table = pd.DataFrame(Hit_List, columns = ['Genome_1', 'Genome_2', 'AAI', 'SD', 'Fragments', 'Total_Fragments_Shortest'])
    # Filter AAI Table by fragment fraction used
    print("Filter by fraction of genome used")
    AAI_Table_Fraction = AAI_Table[AAI_Table['Fragments'] / AAI_Table['Total_Fragments_Shortest'] >= 0.0]
    AAI_Table_Fraction.to_csv(Output, sep="\t", header=True, index=False)
    
    return AAI_Table_Fraction


print("Getting number of genes per genome")
(Number_Genes, Gene_Length) =  Genes_Per_Genome("Protein_DB.faa", "_")

print("Filtering sword by best hit per genome")
BH_per_Genome("02.Test_AAI.sword.search", "02.Test_AAI.sword.search.bh", Gene_Length, "_", "--")

