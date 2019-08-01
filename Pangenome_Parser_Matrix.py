#!/usr/bin/env python

"""
########################################################################
# Author:      Carlos Ruiz, cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# https://github.com/cruizperez/
# Version:    1.0
# Date:      04 May 2019

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

def Pangenome(Matrix, Core = 80, Flexible = 50, Permutations = 100):
    # Read input matrix into a dataframe
    Pangenome_Matrix = pd.read_csv(Matrix, sep="\t")
    Num_Genomes = Pangenome_Matrix.shape[1]
    Record_count = 0
    PG_Dictionary = {}
    # Iterate through permutations.
    for boot in range(1, Num_Genomes+1):
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

    # Change column names to reflec the percentages requested.
    Num_OCs_Pangenome = 'Num_OCs_Pangenome'
    Num_OCs_Flexible = 'Num_OCs_Flexible_' + str(Flexible) + '%'
    Num_OCs_Core = 'Num_OCs_Core_' + str(Core) + '%'

    Pangenome_Table = pd.DataFrame.from_dict(PG_Dictionary, orient = 'index',columns=['Num_Genomes', Num_OCs_Pangenome, Num_OCs_Flexible, Num_OCs_Core])
    Pangenome_Table = Pangenome_Table.sort_values('Num_Genomes')

    return Pangenome_Table


def Pangenome_Plot(Pangenome_Table, Output_Files_Prefix):
    # Import required modules
    from lmfit import minimize, Parameters, fit_report
    from lmfit.models import ExponentialModel, PowerLawModel
    import matplotlib.pyplot as plt


    # Model the curves, we fit an exponential for the pangenome and a powerlaw for the core and flexible

    # For the exponential function.
    #Expo_parameters = Parameters()
    #Expo_parameters.add('amp', value = -1)
    #Expo_parameters.add('decay', value = 0)
    #Expo_parameters.add('offset', value = 0)

    PowerLawModel = PowerLawModel()
    new_param = PowerLawModel.guess(Pangenome_Table.iloc[:,1], x=Pangenome_Table['Num_Genomes'])
    out = PowerLawModel.fit(Pangenome_Table.iloc[:,1], new_param, x=Pangenome_Table['Num_Genomes'])


    ExponentialModel = ExponentialModel()
    new_param2 = ExponentialModel.guess(Pangenome_Table.iloc[:,2], x=Pangenome_Table['Num_Genomes'])
    out2 = ExponentialModel.fit(Pangenome_Table.iloc[:,2], new_param2, x=Pangenome_Table['Num_Genomes'])

    #print(out.fit_report(min_correl=0.25))
    # Define the exponential function
    #def exponenial_func(x, amp, decay, offset):
    #    return amp * np.exp(-decay*x) + offset
    #
    #print(out.plot)
    #print(dir(out))
    ## Define a function to calculate the residual distances for least squares
    #def get_residuals_exponential(params, x, data):
    #    amp = params['amp'].value
    #    decay = params['decay'].value
    #    offset = params['offset'].value
    #
    #    model = exponenial_func(x, amp, decay, offset)
    #
    #    return data - model
    #
    #Output = minimize(get_residuals_exponential, Expo_parameters, args=(Pangenome_Table['Num_Genomes'], Pangenome_Table['Num_OCs_Pangenome']))
    #print(fit_report(out))
    #
    #xx = np.linspace(Pangenome_Table['Num_Genomes'], 200, 500)
    plt.figure(figsize=(10, 10), dpi=300, facecolor='w', edgecolor='k')
    plt.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table['Num_OCs_Pangenome'], alpha = 0.08)
    plt.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table.iloc[:,2], alpha = 0.08, color = 'r')
    plt.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table.iloc[:,3], alpha = 0.08, color = 'g')
    plt.plot(Pangenome_Table['Num_Genomes'], out.best_fit, 'r-')
    plt.plot(Pangenome_Table['Num_Genomes'], out2.best_fit, 'r-')
    plt.savefig(Output_Files_Prefix + '.pdf', dpi=300)
#plt.plot(Pangenome_Table['Num_Genomes'],yy, "k-")

#%%
# Create the parameters of the function.
#from sklearn.metrics import mean_squared_error
#
#def exponenial_func(x, a, b, c):
#    return a*np.exp(-b*x)+c
#
#expon = exponenial_func
#
#def func_powerlaw(x, m, c, c0):
#    return c0 + x**m * c
#
#target_func = func_powerlaw
#
#X = table[0]
#Y = np.array(table[1])
#print(np.shape(Y))
#popt, pcov = curve_fit(target_func, X, table[1], maxfev=2000, p0 = np.asarray([0,10,0]))
#xx = np.linspace(table[0].min(), 200, 100)
#yy = func_powerlaw(X, *popt)
#print("YEAH")
#print(np.shape(yy))
#print(mean_squared_error(Y, yy))
#
#popt2, pcov2 = curve_fit(expon, X, Y, maxfev=2000, p0 = np.asarray([0,10,0]))
#xx2 = np.linspace(table[0].min(), table[0].max(), 100)
#yy2 = exponenial_func(X, *popt2)
#print(mean_squared_error(table[1], yy2))
#
##popt2, pcov2 = curve_fit(expon, table[0], table[2], p0 = np.asarray([1,10,0]), maxfev=5000)
##xx2 = np.linspace(table[0].min(), 150, 100)
##yy2 = exponenial_func(xx2, *popt2)
##
##popt3, pcov3 = curve_fit(target_func, table[0], table[3], p0 = np.asarray([1,10,0]), maxfev=5000)
##xx3 = np.linspace(table[0].min(), 150, 100)
##yy3 = func_powerlaw(xx3, *popt3)
#
##Model = np.polyfit(table[0], table[1], 2)
##plt.subplot(3, 1, 1)
#plt.scatter(table[0], table[1], alpha = 0.01)
##plt.subplot(3, 1, 2)
#plt.scatter(table[0], table[2], color = 'r', alpha = 0.01)
##plt.subplot(3, 1, 3)
#plt.scatter(table[0], table[3], color = 'y', alpha = 0.01)
##plt.show()
##plt.plot(X,yy, "k-")
##plt.plot(xx2,yy2, "b-")
##plt.plot(xx3,yy3, "r-")
#plt.savefig("Pangenome.png")

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
    Pangenome_Table = Pangenome(Input_Matrix, Core, Flexible, Permutations)
    Pangenome_Table.to_csv(Output_Files_Prefix + '.tsv', sep="\t", index=False)

    if Build_Plot == True:
        Pangenome_Plot(Pangenome_Table, Output_Files_Prefix)



if __name__ == "__main__":
    main()
