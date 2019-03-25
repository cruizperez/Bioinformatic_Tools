import numpy as np
import pandas as pd
import seaborn as sns
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
AAI_Matrix = pd.read_csv("04.AAI.matrix", sep="\t", index_col=0)
fig = sns.clustermap(AAI_Matrix, metric='euclidean', method='average', figsize= (20,20))
fig.savefig("output.pdf")
