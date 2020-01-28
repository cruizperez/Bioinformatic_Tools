#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 Nov 27 2019

# Description: Calculates the average amino acid identity using k-mers
from single copy genes. It is a faster version of the regular AAI (Blast
or Diamond) and the hAAI implemented in MiGA.
########################################################################
"""

#! TODO: COMPLETE SCRIPT

# --- Initialize function ---
def child_initialize(_dictionary, _output):
     global Sequence_Dictionary, Output
     Sequence_Dictionary = _dictionary
     Output = _output

def Get_Alignment_Sequences(Alignment_File):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    Dictionary = {}
    with open(Alignment_File) as Fasta:
        for title, sequence in SimpleFastaParser(Fasta):
            Dictionary[title] = sequence
    return Dictionary

def Calculate_Identity(Sequence_ID):
    with open(Output, 'a') as Identity_File:
        for Title_B, Seq_B in Sequence_Dictionary.items():
            Seq_A = Sequence_Dictionary[Sequence_ID]
            Aln_Len = 0
            Match_Len = 0
            Len_Iter = max(len(Sequence_Dictionary[Sequence_ID]), len(Seq_B))
            for i in range(Len_Iter):
                if Seq_A[i] == "-" and Seq_B[i] == "-":
                    continue
                elif Seq_A[i] == Seq_B[i]:
                    Aln_Len += 1
                    Match_Len += 1
                else:
                    Aln_Len += 1
            Identity_File.write("{}\t{}\t{}\n".format(Sequence_ID,Title_B,Match_Len/Aln_Len))


def main():
    import multiprocessing
    Sequences = Get_Alignment_Sequences("04.Archaea_Bacteria.16S.fa.reduced.aln.trim")
    Sequence_IDs = Sequences.keys()

    try:
        pool = multiprocessing.Pool(20, initializer = child_initialize, initargs = (Sequences, "05.16S_Identities.txt"))
        pool.map(Calculate_Identity, Sequence_IDs)
    finally:
        pool.close()
        pool.join()

main()