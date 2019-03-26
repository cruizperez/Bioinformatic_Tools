import pandas as pd
import numpy
import math

Dataframe = pd.read_csv("16S_Matrix.txt", sep="\t", index_col=0)
Header = list(Dataframe.columns.values)
number = len(Header)

for i in range(0, number):
    for j in range(0, number):
        if math.isnan(Dataframe.iloc[i,j]):
            continue
        else:
            Dataframe.iloc[j,i] = Dataframe.iloc[i,j]
