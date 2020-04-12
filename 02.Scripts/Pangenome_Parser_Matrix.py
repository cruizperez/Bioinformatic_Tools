#!/usr/bin/env python

"""
########################################################################
# Author:      Carlos Ruiz, cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# https://github.com/cruizperez/
# Version:    2.0
# Date:      07 Nov 2019

# Description: Parses a binary matrix representing the Pangenome
 of genomes with Orthologous Clusters (OCs) in the rows and genomes
 in columns.

########################################################################
"""

################################################################################
"""---1.0 Import Global Modules---"""
import pandas as pd
import numpy as np

################################################################################
"""---2.0 Define Functions---"""

def Pangenome(Matrix, Core = 80, Flexible = 50, Permutations = 100, Plot=False):
    # Read input matrix into a dataframe
    Pangenome_Matrix = pd.read_csv(Matrix, sep="\t", header=0, index_col=0)
    Num_Genomes = Pangenome_Matrix.shape[1]
    Record_count = 0
    PG_Dictionary = {}
    Statistics_Dictionary = {}
    # Iterate through permutations.
    for boot in range(1, Permutations+1):
        # Sample from 1 to number of genomes provided.
        for step in range(1, Num_Genomes):
            # Randomly subsample the original matrix.
            Subsampled_PG = Pangenome_Matrix.sample(n=step, replace=False, axis=1)
            Sub_PG_Array = Subsampled_PG.to_numpy()
            #Sum the number of genomes in which the Orthologous Cluster is present.
            OCs_Genomes = np.sum(Sub_PG_Array, axis = 1)
            #For the pangenome count as present if OC is present in at least 1 genome.
            Pangenome_OC = (OCs_Genomes > 0).astype(int)
            Num_OCs_PG = np.sum(Pangenome_OC)
            # For the flexible genome, count as present if the OC sum is larger than provided percentage
            Flexible_OC = (((OCs_Genomes*100)/step) > Flexible).astype(int)
            Num_OCs_FG = np.sum(Flexible_OC)
            # For the core genome, count as present if the OC sum is larger than provided percentage
            Core_OC = (((OCs_Genomes*100)/step) > Core).astype(int)
            Num_OCs_CG = np.sum(Core_OC)
            # Add all the values into a list to add in a dictionary.
            lista = [step,Num_OCs_PG,Num_OCs_FG,Num_OCs_CG]
            PG_Dictionary[Record_count] = lista
            Record_count += 1
            # Add values 

    # Change column names to reflec the percentages requested.
    Num_OCs_Pangenome = 'Num_OCs_Pangenome'
    Num_OCs_Flexible = 'Num_OCs_Flexible_' + str(Flexible) + '%'
    Num_OCs_Core = 'Num_OCs_Core_' + str(Core) + '%'

    Pangenome_Table = pd.DataFrame.from_dict(PG_Dictionary, orient = 'index',columns=['Num_Genomes', Num_OCs_Pangenome, Num_OCs_Flexible, Num_OCs_Core])
    Pangenome_Table = Pangenome_Table.sort_values('Num_Genomes')

    if Plot == False:
        import lmfit
        from lmfit.models import ExponentialModel, PowerLawModel, ExpressionModel

        # Model the Pangenome using a Powerlaw function  Ps = κn^γ
        print('\nModelling the Pangenome using a PowerLaw function K*N^\u03B3')
        PowerLawModel = PowerLawModel()
        PowerLaw_Parameters_Guess = PowerLawModel.guess(Pangenome_Table.iloc[:,1], x=Pangenome_Table['Num_Genomes'])
        PowerLaw_Fit = PowerLawModel.fit(Pangenome_Table.iloc[:,1], PowerLaw_Parameters_Guess, x=Pangenome_Table['Num_Genomes'])
        K = float(PowerLaw_Fit.best_values['amplitude'])
        Gamma = float(PowerLaw_Fit.best_values['exponent'])
        print("The equation describing the Pangenome is: {} * N ^ {}".format(round(K,3),round(Gamma,3)))
        if Gamma < 0:
            print("The Pangenome of these organisms is closed with \u03B3 < 0.")
            print("See 10.3389/fmicb.2018.00577 for reference.")
        elif Gamma <= 1:
            print("The Pangenome of these organisms is open with 0 < \u03B3 < 1.")
            print("See 10.3389/fmicb.2018.00577 for reference.")

        # ------------------------

        # Model the Core Genome using an Exponential Decay function  Fc = Kc*exp(-N/τc) + Ω
        print('\nModelling the Core and Flexible Genome using an Exponential Decay function K*exp-(N/\u03C4) + \u03A9')
        
        Custom_Exponential = ExpressionModel('A * exp(-x/tau) + omega', independent_vars=['x'])
        Expression_Parameters = Custom_Exponential.make_params()
        # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
        Expression_Parameters.add_many(('A', 5, True, 0, None, None, 0.1),
                                        ('tau', 5, True, 0, None, None, 0.1),
                                        ('omega', 5, True, 0, None, None, 0.1))
        Core_Genome_Fit = Custom_Exponential.fit(Pangenome_Table.iloc[:,3], Expression_Parameters, x=Pangenome_Table['Num_Genomes'])
        Flexible_Genome_Fit = Custom_Exponential.fit(Pangenome_Table.iloc[:,2], Expression_Parameters, x=Pangenome_Table['Num_Genomes'])

        print("The equation describing the Core and Flexible genomes is: F = K*exp-(N/\u03C4) + \u03A9")
        print("Where \u03A9 represents the approximate number of genes where the estimate reaches a plateau")
        print("Core genome {}% \u03A9: {}".format(Core, round(Core_Genome_Fit.best_values['omega'], 2)))
        print("Flexible genome {}% \u03A9: {}".format(Flexible, round(Flexible_Genome_Fit.best_values['omega'], 2)))

    return Pangenome_Table

# ------------------------------------------------

def Pangenome_Plot(Pangenome_Table, Output_Files_Prefix, Core=80, Flexible=50):
    # Import required modules
    import matplotlib.pyplot as plt
    import lmfit
    from lmfit.models import ExponentialModel, PowerLawModel, ExpressionModel

    # Model the Pangenome using a Powerlaw function  Ps = κn^γ
    print('\nModelling the Pangenome using a PowerLaw function K*N^\u03B3')
    PowerLawModel = PowerLawModel()
    PowerLaw_Parameters_Guess = PowerLawModel.guess(Pangenome_Table.iloc[:,1], x=Pangenome_Table['Num_Genomes'])
    PowerLaw_Fit = PowerLawModel.fit(Pangenome_Table.iloc[:,1], PowerLaw_Parameters_Guess, x=Pangenome_Table['Num_Genomes'])
    K = float(PowerLaw_Fit.best_values['amplitude'])
    Gamma = float(PowerLaw_Fit.best_values['exponent'])
    print("The equation describing the Pangenome is: {} * N ^ {}".format(round(K,3),round(Gamma,3)))
    if Gamma < 0:
        PanLabel = "Closed"
        print("The Pangenome of these organisms is closed with \u03B3 < 0.")
        print("See 10.3389/fmicb.2018.00577 for reference.")
    elif Gamma <= 1:
        PanLabel = "Open"
        print("The Pangenome of these organisms is open with 0 < \u03B3 < 1.")
        print("See 10.3389/fmicb.2018.00577 for reference.")
    # ------------------------
    # Model the Core Genome using an Exponential Decay function  Fc = Kc*exp(-N/τc) + Ω
    print('\nModelling the Core and Flexible Genome using an Exponential Decay function K*exp-(N/\u03C4) + \u03A9')
    
    Custom_Exponential = ExpressionModel('A * exp(-x/tau) + omega', independent_vars=['x'])
    Expression_Parameters = Custom_Exponential.make_params()
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    Expression_Parameters.add_many(('A', 5, True, 0, None, None, 0.1),
                                    ('tau', 5, True, 0, None, None, 0.1),
                                    ('omega', 5, True, 0, None, None, 0.1))
    Core_Genome_Fit = Custom_Exponential.fit(Pangenome_Table.iloc[:,3], Expression_Parameters, x=Pangenome_Table['Num_Genomes'])
    Flexible_Genome_Fit = Custom_Exponential.fit(Pangenome_Table.iloc[:,2], Expression_Parameters, x=Pangenome_Table['Num_Genomes'])
    print("The equation describing the Core and Flexible genomes is: F = K*exp-(N/\u03C4) + \u03A9")
    print("Where \u03A9 represents the approximate number of genes where the estimate reaches a plateau")
    print("Core genome {}% \u03A9: {}".format(Core, round(Core_Genome_Fit.best_values['omega'], 2)))
    print("Flexible genome {}% \u03A9: {}".format(Flexible, round(Flexible_Genome_Fit.best_values['omega'], 2)))
    print("See 10.3389/fmicb.2018.00577 for reference.")
    

    # Get statistics
    Num_Genomes = max(Pangenome_Table['Num_Genomes']) + 1
    # Store lists of Mean, Median, Std
    Pangenome = [[],[],[]]
    Core_Genome = [[],[],[]]
    Flexible_Genome = [[],[],[]]

    for genome in range(1,Num_Genomes):
        Temp_Table = Pangenome_Table[Pangenome_Table.Num_Genomes == genome]
        Pangenome[0].append(np.mean(Temp_Table.iloc[:,1]))
        Pangenome[1].append(np.median(Temp_Table.iloc[:,1]))
        Pangenome[2].append(np.std(Temp_Table.iloc[:,1]))
        Core_Genome[0].append(np.mean(Temp_Table.iloc[:,3]))
        Core_Genome[1].append(np.median(Temp_Table.iloc[:,3]))
        Core_Genome[2].append(np.std(Temp_Table.iloc[:,3]))
        Flexible_Genome[0].append(np.mean(Temp_Table.iloc[:,2]))
        Flexible_Genome[1].append(np.median(Temp_Table.iloc[:,2]))
        Flexible_Genome[2].append(np.std(Temp_Table.iloc[:,2]))


    Figure, Axis = plt.subplots(1,1, figsize=(15,7), constrained_layout=True, dpi=300)
    # Plot scatters
    Axis.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table['Num_OCs_Pangenome'], color="#393373", alpha = 0.2, label="Pangenome", s=70)
    Axis.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table.iloc[:,3], alpha = 0.2, color = '#A82324', label="Core Genome {}%".format(Core), s=70)
    Axis.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table.iloc[:,2], alpha = 0.2, color = '#688170', label="Flexible Genome {}%".format(Flexible), s=70)
    # Plot statistics
    # Axis.scatter(range(1,Num_Genomes), Pangenome[0], color="grey", s=40, marker="s", label="Mean")
    Axis.plot(range(1,Num_Genomes), np.asarray(Pangenome[0]) - np.asarray(Pangenome[2]), linewidth=1, color="black", ls=(0, (1, 1)), label="Standard Dev.")
    Axis.plot(range(1,Num_Genomes), np.asarray(Pangenome[0]) + np.asarray(Pangenome[2]), linewidth=1, color="black", ls=(0, (1, 1)))
    # Axis.scatter(range(1,Num_Genomes), Core_Genome[0], color="grey", s=40, marker="s")
    Axis.plot(range(1,Num_Genomes), np.asarray(Core_Genome[0]) - np.asarray(Core_Genome[2]), linewidth=1, color="black", ls=(0, (1, 1)))
    Axis.plot(range(1,Num_Genomes), np.asarray(Core_Genome[0]) + np.asarray(Core_Genome[2]), linewidth=1, color="black", ls=(0, (1, 1)))
    # Axis.scatter(range(1,Num_Genomes), Flexible_Genome[0], color="grey", s=40, marker="s")
    Axis.plot(range(1,Num_Genomes), np.asarray(Flexible_Genome[0]) - np.asarray(Flexible_Genome[2]), linewidth=1, color="black", ls=(0, (1, 1)))
    Axis.plot(range(1,Num_Genomes), np.asarray(Flexible_Genome[0]) + np.asarray(Flexible_Genome[2]), linewidth=1, color="black", ls=(0, (1, 1)))
    # Plot best fits
    Axis.plot(Pangenome_Table['Num_Genomes'], PowerLaw_Fit.best_fit, '-', linewidth=2, label="\u03B3: {} ({} Pangenome)".format(round(Gamma,2),PanLabel), color="#BFAA2E")
    Axis.plot(Pangenome_Table['Num_Genomes'], Core_Genome_Fit.best_fit, '-', linewidth=2, color="#015C1B", label="Core genes - \u03A9: {}".format(round(Core_Genome_Fit.best_values['omega'])))
    Axis.plot(Pangenome_Table['Num_Genomes'], Flexible_Genome_Fit.best_fit, '-', linewidth=2, color="#824D55", label="Common genes - \u03A9: {}".format(round(Flexible_Genome_Fit.best_values['omega'])))
    Legend = Axis.legend(fontsize=12)
    for text in Legend.get_texts():
        text.set_color("grey")
    Axis.tick_params(axis='both', which='major', labelsize=10)
    Axis.set_xlabel("Number of Genomes", fontsize=15)
    Axis.set_ylabel("Number of Gene Clusters", fontsize=15)
    Axis.set_xlim(0,Num_Genomes)
    Axis.set_ylim(min(Pangenome_Table.iloc[:,3])-500,max(Pangenome_Table.iloc[:,1])+500)
    Figure.suptitle("Core vs Pangenome Collector Curves", color="#802B18", fontsize=35)
    Figure.savefig(Output_Files_Prefix + ".png")

# ----------------------------------------------


"""---3.0 Main Function---"""

def main():
    from sys import argv
    import argparse
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Parses a binary matrix representing the Pangenome
                                                    of genomes with Orthologous Clusters (OCs) in the rows and genomes
                                                    in columns. It produces a tab-delimited file with the Pan, Core and Flexible genomes and optionally a plot'''
                                    'Global mandatory parameters: [Input Binary Matrix] [Output Prefix]\n'
                                    'Optional Database Parameters: See ' + argv[0] + ' -h')
    parser.add_argument("-i", "--inputFile", dest='Input_Matrix', required=True, help="Input FastA file")
    parser.add_argument("-o", "--outputPrefix", dest='Output_Files_Prefix', required=False, default = 'Pangenome', help="Prefix for the output files.")
    parser.add_argument("--core", dest='Core', required=False, type = int, default = 80, help="Percentage of genomes with an OC to be considered core, by default 80")
    parser.add_argument("--flex", dest='Flexible', required=False, type = int, default = 50, help="Percentage of genomes with an OC to be considered part of the flexible genome, by default 50")
    parser.add_argument("--plot", dest='Build_Plot', required=False, action='store_true', help="True/False. Build the plot with the Pan, Core and Flexible genomes, by default false")
    parser.add_argument("--permutation", dest='Permutations', required=False, type = int, default= 100, help="Number of permutations to execute, by default 100")

    args = parser.parse_args()

    # Parse input arguments
    Input_Matrix = args.Input_Matrix
    Output_Files_Prefix = args.Output_Files_Prefix
    Core = args.Core
    Flexible = args.Flexible
    Build_Plot = args.Build_Plot
    Permutations = args.Permutations

    # Run matrix Parser and save matrix
    Pangenome_Table = Pangenome(Input_Matrix, Core, Flexible, Permutations, Build_Plot)
    Pangenome_Table.to_csv(Output_Files_Prefix + '.tsv', sep="\t", index=False)

    if Build_Plot == True:
        Pangenome_Plot(Pangenome_Table, Output_Files_Prefix, Core, Flexible)



if __name__ == "__main__":
    main()
