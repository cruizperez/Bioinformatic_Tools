#!/usr/bin/env python

"""----------------------------- 1.0 Define Functions -----------------------------"""

def FastA_Filter_List(FastaFile, Output, List, Reverse=False):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    Seq_ID_list = []
    if type(List) == list:
        Seq_ID_list = List
        Records = len(Seq_ID_list)
    else:
        import pandas as pd
        List_DF = pd.read_csv(List, delim_whitespace=True)
        Seq_ID_list = df[0].tolist()
    Records = len(Seq_ID_list)
    with open(FastaFile) as Fasta_in, open(Output, 'w') as Fasta_out:
        if Reverse == True:
            print("Excluding " + str(Records) + " records from output")
            for title, seq in SimpleFastaParser(Fasta_in):
                if title.split()[0] not in Seq_ID_list:
                    Fasta_out.write(">%s\n%s\n" % (title, seq))
        else:
            print("Retrieving " + str(Records) + " records from input")
            for title, seq in SimpleFastaParser(Fasta_in):
                if len(Seq_ID_list) < 1:
                    break
                elif title.split()[0] in Seq_ID_list:
                    Fasta_out.write(">%s\n%s\n" % (title, seq))
                    Seq_ID_list.remove(title.split()[0])
    Fasta_out.close()


def main():
    import argparse
    from sys import argv
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description='''Filter a FastA file based ona provided list of IDs (file or input).\n'''
    '''It can exclude or retrieve the sequences using the --reverse flag\n'''
    '''Usage: ''' + argv[0] + ''' -f [FastA File] -o [Output File] -l [File with IDs] OR -i [ID list]\n'''
    '''Global mandatory parameters: [FastA_File] [Output_File] [ID List File]\n'''
    'Optional Database Parameters: See ' + argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='Fasta_File', action='store', required=True, help='FastA file to filter')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output FastA file with retrieved sequences')
    parser.add_argument('-l', '--list', dest='ID_File', action='store', required=False, help='File with list of IDs to filter.')
    parser.add_argument('-i', '--id', dest='ID_List', action='store', required=False, help='Comma-separated IDs to filter: ID1,ID2,ID3')
    parser.add_argument('--reverse', action='store_true', help='Exclude the sequences in the list file. By default False, i.e. retrieves those in the list')
    args = parser.parse_args()

    Fasta_File = args.Fasta_File
    Output_File = args.Output_File
    ID_File = args.ID_File
    ID_List = args.ID_List
    Reverse = args.reverse

    if ID_List == None:
        FastA_Filter_List(Fasta_File, Output_File, ID_File, Reverse)
    elif ID_File == None:
        ID_List = ID_List.split(",")
        FastA_Filter_List(Fasta_File, Output_File, ID_List, Reverse)
    else:
        raise ValueError("Did you provide the IDs to filter?")

if __name__ == "__main__":
    main()
