#!/usr/bin/env python

"""
########################################################################
# Author:      Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:    1.0
# Date:      04 May 2019

# Description: This script creates a MySQL database from the NCBI bacterial,
# archaeal and viral RefSeq databases for annotation purposes.

########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import os, sys, argparse
import mysql.connector
from Bio import SeqIO



#################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def Pangenome(Matrix, Core = 80, Flexible = 50, Bootstraps = 100):
    # Read input matrix into a dataframe
    Pangenome_Matrix = pd.read_csv(Matrix, sep="\t")
    Num_Genomes = Pangenome_Matrix.shape[1]
    Record_count = 0
    PG_Dictionary = {}
    # Iterate through bootstraps.
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


Pangenome_Table = Pangenome("02.Pangenome_Analysis.txt", 70, 20, 100)
Pangenome_Table.to_csv("pangenome.txt", sep="\t", index=False)

#%%

from lmfit import minimize, Parameters, fit_report
from lmfit.models import ExponentialModel, PowerLawModel
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
plt.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table['Num_OCs_Pangenome'], alpha = 0.01)
plt.scatter(Pangenome_Table['Num_Genomes'], Pangenome_Table.iloc[:,2], alpha = 0.01, color = 'r')
plt.plot(Pangenome_Table['Num_Genomes'], out.best_fit, 'r-')
plt.plot(Pangenome_Table['Num_Genomes'], out2.best_fit, 'r-')
#plt.plot(Pangenome_Table['Num_Genomes'],yy, "k-")

#%%
# Create the parameters of the function.
from sklearn.metrics import mean_squared_error

def exponenial_func(x, a, b, c):
    return a*np.exp(-b*x)+c

expon = exponenial_func

def func_powerlaw(x, m, c, c0):
    return c0 + x**m * c

target_func = func_powerlaw

X = table[0]
Y = np.array(table[1])
print(np.shape(Y))
popt, pcov = curve_fit(target_func, X, table[1], maxfev=2000, p0 = np.asarray([0,10,0]))
xx = np.linspace(table[0].min(), 200, 100)
yy = func_powerlaw(X, *popt)
print("YEAH")
print(np.shape(yy))
print(mean_squared_error(Y, yy))

popt2, pcov2 = curve_fit(expon, X, Y, maxfev=2000, p0 = np.asarray([0,10,0]))
xx2 = np.linspace(table[0].min(), table[0].max(), 100)
yy2 = exponenial_func(X, *popt2)
print(mean_squared_error(table[1], yy2))

#popt2, pcov2 = curve_fit(expon, table[0], table[2], p0 = np.asarray([1,10,0]), maxfev=5000)
#xx2 = np.linspace(table[0].min(), 150, 100)
#yy2 = exponenial_func(xx2, *popt2)
#
#popt3, pcov3 = curve_fit(target_func, table[0], table[3], p0 = np.asarray([1,10,0]), maxfev=5000)
#xx3 = np.linspace(table[0].min(), 150, 100)
#yy3 = func_powerlaw(xx3, *popt3)

#Model = np.polyfit(table[0], table[1], 2)
#plt.subplot(3, 1, 1)
plt.scatter(table[0], table[1], alpha = 0.01)
#plt.subplot(3, 1, 2)
plt.scatter(table[0], table[2], color = 'r', alpha = 0.01)
#plt.subplot(3, 1, 3)
plt.scatter(table[0], table[3], color = 'y', alpha = 0.01)
#plt.show()
#plt.plot(X,yy, "k-")
#plt.plot(xx2,yy2, "b-")
#plt.plot(xx3,yy3, "r-")
plt.savefig("Pangenome.png")
