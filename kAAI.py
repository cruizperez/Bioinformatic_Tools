#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 Oct 04 2019

# Description: This script parses a CheckM "bin_stats_ext.tsv" file from
the analyze workflow along with a list of SCGs and returns a matrix with 
each genome and SCG found in it.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""
def run_prodigal(InputFile):
    import subprocess
    from pathlib import Path

    FilePath = Path(InputFile)
    Prefix = Path(FilePath.stem)
    Folder = FilePath.parent
    Output = Folder / Prefix.with_suffix('.faa')
    Temp_Output = Folder / Prefix.with_suffix('.temp')
    subprocess.call(["prodigal", "-i", str(FilePath), "-a", str(Output), "-p", "meta", "-q", "-o", str(Temp_Output)])
    Temp_Output.unlink()
    return Output

def run_hmmsearch(InputFile):
    import subprocess
    from pathlib import Path

    FilePath = Path(InputFile)
    Prefix = Path(FilePath.stem)
    Folder = FilePath.parent
    Output = Folder / Prefix.with_suffix('.hmm')
    Temp_Output = Folder / Prefix.with_suffix('.temp')
    subprocess.call(["hmmsearch", "--tblout", str(Output), "-o", str(Temp_Output), "--cut_ga", "--cpu", "1", "/mnt/c/Users/Cruiz/OneDrive/Escritorio/kAAI/Multiprocess_Test/00.HMM_Models/genes.hmm", str(FilePath)])
    Temp_Output.unlink()
    return Output

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

def read_kmers_from_file(filename, positive_hits, ksize):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    all_kmers = []
    with open(filename) as Fasta_in:
        for title, sequence in SimpleFastaParser(Fasta_in):
            if title.split()[0] in positive_hits:
                kmers = build_kmers(sequence, ksize)
                all_kmers += kmers
    return all_kmers

def kAAI_Parser(ID, Kmer_Dictionary, OutFile):
    with open(OutFile, 'a') as OutFile:
        # for ID in ID_List:
        for key2, value2 in Kmer_Dictionary.items():
            intersection = len(set(Kmer_Dictionary[ID]).intersection(set(value2)))
            shorter = min(len(list(set(Kmer_Dictionary[ID]))), len(list(set(value2))))
            fraction = round(intersection/shorter, 3)
            OutFile.write("{}\t{}\t{}\t{}\t{}\n".format(ID, key2, intersection, shorter, fraction))

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse
    from sys import argv
    from sys import exit
    from pathlib import Path
    import subprocess
    import multiprocessing
    from functools import partial

    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a CheckM "bin_stats_ext.tsv" file from\n'''
                        '''the analyze workflow along with a list of SCGs and returns a matrix with\n'''
                        '''each genome and SCG found in it.\n'''
            '''Usage: ''' + argv[0] + ''' -i [CheckM Output] -l [SCG list] -o [Gene Copy Matrix]\n'''
            '''Global mandatory parameters: -i [CheckM Output] -l [SCG list] -o [Gene Copy Matrix]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    parser.add_argument('-g', '--genomes', dest='Genome_List', action='store', nargs='+', required=False, help='List of input genomes')
    parser.add_argument('-p', '--proteins', dest='Protein_Files', action='store', nargs='+', required=False, help='List of input protein files')
    parser.add_argument('-o', '--output', dest='Output', action='store', required=True, help='Output File')
    parser.add_argument('-t', '--threads', dest='Threads', action='store', default=1, type=int, required=False, help='Number of threads to use, by default 1')
    args = parser.parse_args()

    Genome_List = args.Genome_List
    Protein_Files = args.Protein_Files
    Output = args.Output
    Threads = args.Threads

    if Genome_List != None:
        print("Predicting proteins...")
        Protein_Files = []
        HMM_Search_Files = []
        for Genome in Genome_List:
            Proteins = run_prodigal(Genome)
            Protein_Files.append(Proteins)
        print("Searching HMM models...")
        for Proteins in Protein_Files:
            HMM_Result = run_hmmsearch(Proteins)
            HMM_Search_Files.append(HMM_Result)
    elif Protein_Files != None:
        print("Searching HMM models...")
        number_processes = Threads
        pool = multiprocessing.Pool(number_processes)
        HMM_Search_Files = pool.map(run_hmmsearch, Protein_Files)
        pool.close()
        pool.join()
    else:
        exit('No input provided, please provide either genomes "-g" or protein files "-p"')

    Kmer_Dic = {}
    for SCG_file in Protein_Files:
        Protein_Path = Path(SCG_file)
        Prefix = Path(Protein_Path.stem)
        Name = Protein_Path.name
        Folder = Protein_Path.parent
        HMM_File = Folder / Prefix.with_suffix('.hmm')
        Positive_Matches = []
        with open(HMM_File, 'r') as HMM_Input:
            for line in HMM_Input:
                if line.startswith("#"):
                    continue
                else:
                    Positive_Matches.append(line.strip().split()[0])
        HMM_File.unlink()
        kmers = read_kmers_from_file(SCG_file, Positive_Matches, 8)
        Kmer_Dic[Name] = kmers
    
    ID_List = Kmer_Dic.keys()
    try:
        pool = pool = multiprocessing.Pool(number_processes)
        pool.map(partial(kAAI_Parser, Kmer_Dictionary=Kmer_Dic, OutFile=Output), ID_List)
    finally:
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()