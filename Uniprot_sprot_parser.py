#!/usr/bin/env python

from collections import defaultdict
import pandas as pd

dictionary = defaultdict(list)
identification = None

with open("uniprot_trembl.dat") as file:
    for line in file:
        line = line.strip()
        line = line.split("   ")
        if len(line) <= 1:
            continue
        elif line[0] == "ID":
            identification = line[1]
            dictionary[identification] = [[] for k in range(5)]
        elif line[0] == "AC":
            accession = line[1]
            dictionary[identification][0].extend([accession])
        elif line[0] == "DE" and "RecName" in line[1]:
            recname = line[1].split("=")[1]
            dictionary[identification][1].extend([recname])
        elif line[0] == "DR" and "GO:" in line[1]:
            if "C:" in line[1]:
                compartment = line[1].split("C:")[1]
                dictionary[identification][2].extend([compartment])
            elif "F:" in line[1]:
                function = line[1].split("F:")[1]
                dictionary[identification][3].extend([function])
            elif "P:" in line[1]:
                process = line[1].split("P:")[1]
                dictionary[identification][4].extend([process])


for key in dictionary:
    dictionary[key][0] = [' '.join(dictionary[key][0])]
    dictionary[key][1] = [' '.join(dictionary[key][1])]
    dictionary[key][2] = [' '.join(dictionary[key][2])]
    dictionary[key][3] = [' '.join(dictionary[key][3])]
    dictionary[key][4] = [' '.join(dictionary[key][4])]

df=pd.DataFrame.from_dict(dictionary,orient='index')
df.columns = ["Accesion", "Gene", "Compartment", "Function", "Process"]
for column in df:
    df[column] = df[column].str[0]
df.to_csv("Salida.txt", sep="\t")
