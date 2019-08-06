#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 August 5 2019

# Description: This script simulates a microbial community based on FastA file and a complexity
based on a power law model, teh choices are low, medium and high complexities.
Parameters based on Johnson, S., Trost, B., Long, J. R., Pittet, V., & Kusalik, A. (2014).
A better sequence-read simulator program for metagenomics. BMC bioinformatics, 15 Suppl 9(Suppl 9), S14.
doi:10.1186/1471-2105-15-S9-S14
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Community_Simulation(Input_File, Output_Name, Complexity, Min_Species, Max_Species, Constant, Exponent):
    from random import randint
    import pandas as pd
    from random import sample
    # Create empty dataframe
    Community_Abundances = pd.DataFrame(columns=["Genome", Output_Name])
    # Maximum number of species to simulate
    Max_Species = Max_Species
    Min_Species = Min_Species
    # Read FastA File
    Candidate_Genomes = []
    Input_Genomes = 0
    with open(Input_File) as FastA:
        Input_Genomes = 0
        for line in FastA:
            if ">" in line:
                Input_Genomes += 1
                Candidate_Genomes.append(line.split()[0].replace(">",""))
    # Check maximum number of species wanted vs genomes provided
    if Input_Genomes < Max_Species:
        print("Maximum number of species exceeds the number of genomes\n" +
              "in your input FastA file.\n" +
              "Maximum number of species is now {}.".format(Input_Genomes))
        Max_Species = Input_Genomes
    if Input_Genomes < Min_Species:
        print("Minimum number of species exceeds the number of genomes\n" +
              "in your input FastA file.\n" +
              "Minimum number of species is now 5.")
        Min_Species = 5
    # Select total number of members to simulate
    Total_Members = randint(Min_Species, Max_Species)
    # List selected genomes for simulation
    Selected_Genomes = sample(Candidate_Genomes, Total_Members)

    print("Simulating a microbial community with {} members and a {} complexity"
          .format(Total_Members, Complexity.lower()))

    # Determine the total abundance and the relative abundances.
    Abundance = 0
    Abundance_Total = 0
    Species_Num = 1

    for i in range(1, len(Selected_Genomes)+1):
        Abundance_Total += (Constant * (i**Exponent))

    for i in range(1, len(Selected_Genomes)+1):
        Abundance = (Constant * (Species_Num ** Exponent))
        Abundance = Abundance / Abundance_Total
        Community_Abundances.at[i, "Genome"] = Selected_Genomes[i-1]
        Community_Abundances.at[i, Output_Name] = Abundance
        Species_Num += 1
    #Community_Abundances.set_index('Genome')
    return Community_Abundances

def Add_Community_Plot(Dataframe, Iteration, Complexity, Axes):
    Complex = str(Complexity)
    #print(Dataframe)
    #Serie = Dataframe.iloc[:0]
    #print(Serie)
    Index = "Community_{}".format(Iteration)
    Serie = Dataframe.sort_values(Index, ascending=False)
    Serie = Serie[Serie > 0]
    Serie = Serie.reset_index(drop=True)
    Genomes = len(Serie)
    Axes[Iteration-1].scatter(Serie.index, Serie)
    Axes[Iteration-1].plot(Serie.index, Serie, '--', color='orange',
                label="Number of Genomes: {}\nComplexity: {}".format(str(Genomes), Complex))
    Axes[Iteration-1].title.set_text("Community {}".format(str(Iteration)))
    Axes[Iteration-1].legend(loc='center left', bbox_to_anchor=(0.65, 0.9))
    Axes[Iteration-1].set_ylabel('Relative Abundance')
    Axes[Iteration-1].set_xlabel('Genome Rank')



################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys, os
    from FastA_Filter_List import FastA_Filter_List
    import pandas as pd
    from random import sample
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script simulates a microbial community following a power law distribution.\n'''
            '''All you need to provide is a FastA file, the number of species you want and the complexity\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [FastA File] -o [Output Prefix] -c [Complexity] --min [Minimum # species]
            --max [Maximum # Species] --iterations [Number of communities to simulate] --plot\n'''
            '''Global mandatory parameters: -i [FastA File] -o [Output Prefix] -c [Complexity]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--inputFastA', dest='Input_FastA', action='store', required=True, help='Input FastA to simulate')
    parser.add_argument('-o', '--outputPrefix', dest='Output_Prefix', action='store', required=True, help='Folder to store outputs')
    parser.add_argument('-c', '--complexity', dest='Complexity', action='store', required=True, help='Community complexity, random, low, medium or high')
    parser.add_argument('--min', dest='Min_Species', action='store', required=False, type=int, help='Minimum number of species', default=5)
    parser.add_argument('--max', dest='Max_Species', action='store', required=False, type=int, help='Maximum number of species', default=1000)
    parser.add_argument('--iterations', dest='Iterations', type=int, action='store', required=False, help='Number of communities to simulate', default=1)
    parser.add_argument('--plot', dest='Plot', action='store_true', required=False, help='Plot community rank plot, disabled by default')
    args = parser.parse_args()

    Input_FastA = args.Input_FastA
    Output_Prefix = args.Output_Prefix
    Complexity = args.Complexity
    Min_Species = args.Min_Species
    Max_Species = args.Max_Species
    Iterations = args.Iterations
    Plot = args.Plot

    # ---------------------------------
    # Create output folder if doesn't exist
    if os.path.exists(Output_Prefix) == False:
        os.mkdir(Output_Prefix)
    # Initialize empty dataframe
    Final_Table = pd.DataFrame()

    # Initialize plot if activated
    if Plot == True:
        import matplotlib.pyplot as plt
        Num_Plots = Iterations
        Height = 5 * Num_Plots
        Figure, Axes = plt.subplots(Num_Plots, 1, figsize=(10,Height), sharex=True, sharey = True,
            constrained_layout=True)
        Figure.suptitle('Communities Simulated', fontsize=16)

    for i in range(1,Iterations+1):
        print("\nIteration {}".format(str(i)))
        # Select complexity value
        if Complexity.lower() == "random":
            Complexity_Community = sample(["low", "medium", "high"], 1)[0]
            if Complexity_Community.lower() == "low":
                Complexity_Community = "low"
                Constant = 31.4034355
                Exponent = -1.2868576
            elif Complexity_Community.lower() == "medium":
                Complexity_Community = "medium"
                Constant = 21.2301963
                Exponent = -1.0582662
            elif Complexity_Community.lower() == "high":
                Complexity_Community = "high"
                Constant = 2.08348893
                Exponent = -0.233125
            else:
                raise ValueError("Please specify a level of complexity, i.e. random, low, medium or high")
        else:
            if Complexity.lower() == "low":
                Complexity_Community = "low"
                Constant = 31.4034355
                Exponent = -1.2868576
            elif Complexity.lower() == "medium":
                Complexity_Community = "medium"
                Constant = 21.2301963
                Exponent = -1.0582662
            elif Complexity.lower() == "high":
                Complexity_Community = "high"
                Constant = 2.08348893
                Exponent = -0.233125
            else:
                raise ValueError("Please specify a level of complexity, i.e. random, low, medium or high")

        # Simulate communities
        Community_Name = "Community_" + str(i)
        Abundances = Community_Simulation(Input_FastA, Community_Name, Complexity_Community, Min_Species, Max_Species, Constant, Exponent)
        Sequences_IDs = Abundances['Genome'].tolist()
        Abundances.set_index('Genome', inplace=True)
        if Plot == True:
            Add_Community_Plot(Abundances, i, Complexity_Community, Axes)
        # Filter original fasta with only those members present in the community
        FastA_Out = Output_Prefix + '/' + Community_Name + '.fasta'
        FastA_Filter_List(Input_FastA, FastA_Out, Sequences_IDs)
        if len(Final_Table) == 0:
            Final_Table = Abundances
        else:
            Final_Table = pd.concat([Final_Table, Abundances], axis=1, sort=True)
    # Save final abundance table
    Table_Out = Output_Prefix + '/Community_Abundance.tsv'
    Final_Table = Final_Table.fillna(value=0)
    Final_Table.index.name = 'Genome_ID'
    Final_Table.to_csv(Table_Out, header = True, index = True, sep="\t")

    # Save final figure
    Figure.savefig(Output_Prefix + '/Simulation_Plot.pdf')


if __name__ == "__main__":
    main()
