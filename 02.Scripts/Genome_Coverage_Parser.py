import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


Clades = ["Ia", "Ib", "Ic", "IIa.A", "IIa.B", "IIb", "IIIa", "IIIb", "V"]
Equivalent_Metagenome = {}
Clade_Coverage = pd.DataFrame(columns=["Ia", "Ib", "Ic", "IIa.A", "IIa.B", "IIb", "IIIa", "IIIb", "V"])

with open("01.Metagenome_Genome_Equivalent.txt") as Equivalents:
    for line in Equivalents:
        line = line.strip().split()
        Equivalent_Metagenome[line[0]] = float(line[1])

with open("06.Merged_Coverages.txt") as Coverage:
    for line in Coverage:
        line = line.strip().split()
        if line[0] not in Clade_Coverage.index:
            Clade_Coverage = Clade_Coverage.reindex(Clade_Coverage.index.values.tolist()+[line[0]], fill_value=0)
            for clade in Clades:
                if line[1].startswith(clade):
                    Clade_Coverage.at[line[0], clade] += float(line[2])/Equivalent_Metagenome[line[0]]*100
                else:
                    pass
        else:
            for clade in Clades:
                if line[1].startswith(clade):
                    Clade_Coverage.at[line[0], clade] += float(line[2])/Equivalent_Metagenome[line[0]]*100
                else:
                    pass

#print (Clade_Coverage)

colors = ["#A0A053", "#638960", "#1D4763", "#A02048", "#C78F44", "#BF6459", "#5D676D", "#775E7A", "#6B5856"]
Clade_Coverage.loc[:,["Ia", "Ib", "Ic", "IIa.A", "IIa.B", "IIb", "IIIa", "IIIb", "V"]].plot.barh(stacked=True, figsize=(8,10), color=colors, width=1)

#plt.savefig('books_read.pdf')
