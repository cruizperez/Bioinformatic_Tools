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
# --- Run prodigal ---
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

# --- Run hmmsearch ---
def run_hmmsearch(InputFile):
    import subprocess
    from pathlib import Path

    FilePath = Path(InputFile)
    Prefix = Path(FilePath.stem)
    Folder = FilePath.parent
    Output = Folder / Prefix.with_suffix('.hmm')
    Temp_Output = Folder / Prefix.with_suffix('.temp')
    Script_path = Path(__file__)
    Script_dir = Script_path.parent
    HMM_Model = Script_dir / "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"
    subprocess.call(["hmmsearch", "--tblout", str(Output), "-o", str(Temp_Output), "--cut_ga", "--cpu", "1", str(HMM_Model), str(FilePath)])
    Temp_Output.unlink()
    return Output

# --- Build Kmers ---
def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

# --- Read Kmers from files ---
def read_kmers_from_file(filename, positive_hits, ksize):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    all_kmers = []
    with open(filename) as Fasta_in:
        for title, sequence in SimpleFastaParser(Fasta_in):
            if title.split()[0] in positive_hits:
                kmers = build_kmers(sequence, ksize)
                all_kmers += kmers
    return all_kmers

# --- Read Kmers from files ---
def read_total_kmers_from_file(filename, positive_hits, ksize):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    all_kmers = []
    with open(filename) as Fasta_in:
        for title, sequence in SimpleFastaParser(Fasta_in):
            kmers = build_kmers(sequence, ksize)
            all_kmers += kmers
    return all_kmers

# --- Parse kAAI ---
def kAAI_Parser(ID):
    from pathlib import Path

    FilePath = Path(ID)
    Folder = Path.cwd()
    Output = Folder / FilePath.with_suffix('.aai.temp')
    with open(Output, 'w') as OutFile:
        for key2, value2 in Kmer_Dictionary.items():
            intersection = len(Kmer_Dictionary[ID].intersection(value2))
            shorter = min(len(list(Kmer_Dictionary[ID])), len(list(value2)))
            fraction = round(intersection/shorter, 3)
            OutFile.write("{}\t{}\t{}\t{}\t{}\n".format(ID, key2, intersection, shorter, fraction))
    return Output

# --- Initialize function ---
def child_initialize(_dictionary):
     global Kmer_Dictionary
     Kmer_Dictionary = _dictionary


def Kmer_Parser(SCG_file):
    from pathlib import Path

    Kmer_Dic = {}
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
    Kmer_Dic[Name] = set(kmers)

    return Kmer_Dic

def merge_dicts(Dictionaries):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in Dictionaries:
        result.update(dictionary)
    return result

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
    import datetime
    import shutil

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

    print(datetime.datetime.now()) # Remove after testing
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
        try:
            pool = multiprocessing.Pool(Threads)
            HMM_Search_Files = pool.map(run_hmmsearch, Protein_Files)
        finally:
            pool.close()
            pool.join()
    else:
        exit('No input provided, please provide either genomes "-g" or protein files "-p"')
    print(datetime.datetime.now()) # Remove after testing

    print(datetime.datetime.now()) # Remove after testing
    
    # Parse HMM results
    print("Parsing HMM results...")
    try:
        pool = multiprocessing.Pool(Threads)
        Kmer_Results = pool.map(Kmer_Parser, Protein_Files)
    finally:
        pool.close()
        pool.join()

    Final_Kmer_Dict = merge_dicts(Kmer_Results)
    print(datetime.datetime.now()) # Remove after testing

    # Calculate shared Kmer fraction
    print(datetime.datetime.now()) # Remove after testing
    print("Calculating shared Kmer fraction...")
    ID_List = Final_Kmer_Dict.keys()
    try:
        pool = multiprocessing.Pool(Threads, initializer = child_initialize, initargs = (Final_Kmer_Dict,))
        Fraction_Results = pool.map(kAAI_Parser, ID_List)
        # Fraction_Results = pool.map(partial(kAAI_Parser, Kmer_Dictionary=Final_Kmer_Dict.copy()), ID_List)
    finally:
        pool.close()
        pool.join()
    print(datetime.datetime.now()) # Remove after testing

     # Merge results into a single output
    with open(Output, 'w') as OutFile:
        for file in Fraction_Results:
            with open(file) as Temp:
                shutil.copyfileobj(Temp, OutFile)
            file.unlink()


if __name__ == "__main__":
    main()