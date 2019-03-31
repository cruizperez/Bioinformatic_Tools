def Annotation_Parser(Input_files, Gene_Col, Annotation_Col):
    Annotations = {}
    for file in Input_files:
        with open(file) as Annotation_File:
            for line in Annotation_File:
                line = line.strip().split()
                if line[Gene_Col-1] not in Annotations:
                    try:
                        Annotations[line[Gene_Col-1]] = line[Annotation_Col-1]
                    except:
                        Annotations[line[Gene_Col-1]] = ""
                else:
                    if Annotations[line[Gene_Col-1]] == "":
                        Annotations[line[Gene_Col-1]] = line[Annotation_Col-1]
                    else:
                        pass
    return(Annotations)

Diccion = Annotation_Parser(["13.SAR11_KO_SBH.txt", "14.SAR11_KO_MAPLE.txt"], 1, 3)
Output = open("SALIDA.txt", 'w')
for key, value in Diccion.items():
    Output.write("%s\t%s\n" % (key, value))
Output.close()
