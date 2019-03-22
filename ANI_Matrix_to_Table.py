#!/usr/bin/env python

################################################################################
"""---1.0 Import Modules---"""

import pandas as pd


################################################################################
"""---2.0 Define Functions---"""

def ANI_Parser(ANI_File, Coverage_File, Output_File, ID):
    Coverage_DF = pd.read_csv(Coverage_File, sep="\t", index_col=0)
    Identity_DF = pd.read_csv(ANI_File, sep="\t", index_col=0)

    Genomes = list(Coverage_DF)

    Output_FH = open(Output_File, "w")

    for i in Genomes:
        for j in Genomes:
            if Coverage_DF.loc[i, j] >= 0.2 or Coverage_DF.loc[j, i] >= 0.2 and Identity_DF.loc[i, j] >= ID:
                Output_FH.write("%s\t%s\t%s\n" % (i, j, Identity_DF.loc[i, j]))

    Output_FH.close()

################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(description='''Parses pyANI output returning ANI values that were calculated over more than 20%% of the sequence, if no ID value is provided the default is 95%%'''
                                    'Global mandatory parameters: [ANI_Matrix] [ANI_Aln_Percent_Matrix] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-a', '--ani', dest='ANI_File', action='store', required=True, help='Matrix with ANI values')
    parser.add_argument('-c', '--coverage', dest='Coverage_File', action='store', required=True, help='Matrix with coverage percentage')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output Table file with retrieved ANI pairs and IDs')
    parser.add_argument('-i', '--id', dest='ID', action='store', help='ANI ID Percentage to filter', type=float, default=0.95)
    args = parser.parse_args()

    ANI_File = args.ANI_File
    Coverage_File = args.Coverage_File
    Output_File = args.Output_File
    ID = args.ID

    ANI_Parser(ANI_File, Coverage_File, Output_File, ID):

if __name__ == "__main__":
    main()
