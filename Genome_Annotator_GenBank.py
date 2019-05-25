#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 19 April 2019

# Description: Combines a genome sequence in fastA format, a gff and protein files derived from prodigal,
and optionally a tab-separated annotation file into a single genbank file.
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""
import argparse, sys
from Bio import SeqIO
from BCBio import GFF
from Bio.Alphabet import generic_dna, generic_protein
from os import remove
import datetime

################################################################################
"""---2.0 Define Functions---"""

def GFF_Fasta_Merger(Fasta_Input, GFF_File):
    Genome = SeqIO.to_dict(SeqIO.parse(Fasta_Input, "fasta", generic_dna))
    Genbank_File = GFF.parse(GFF_File, Genome)
    SeqIO.write(Genbank_File, "TEMP_GENBANK.gb", "genbank")

def Protein_2_GenBank(Genbank_Record, Protein_File):
    with open(Protein_File) as proteins:
        for protein in SeqIO.parse(proteins, "fasta"):
            tag = protein.name
            ids = protein.description.split(sep="#")[4].split(sep=";")[0].split(sep="=")[1]
            seq = protein.seq
            for feature in Genbank_Record.features:
                if feature.qualifiers['ID'][0] == ids:
                    feature.qualifiers['locus_tag'] = [tag]
                    feature.qualifiers['translation'] = [str(seq)]
    return Genbank_Record

def Annotation_2_GenBank(Genbank_Record, Annotation_File, NameCol, AnnotCol):
    with open(Annotation_File) as annotations:
        for line in annotations:
            line = line.strip().split(sep="\t")
            for feature in Genbank_Record.features:
                if feature.qualifiers['locus_tag'] == [line[NameCol]]:
                    feature.qualifiers['product'] = [line[AnnotCol]]
    return Genbank_Record


################################################################################
"""---3.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''Combines a genome sequence in fastA format, a gff and protein files derived from prodigal,
    and optionally a tab-separated annotation file into a single genbank file.'''
                                    'Global mandatory parameters: [Genome_FastA] [GFF Gene Prediction] [Protein Fasta]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument("-i", "--inputGenome", dest='Input_Genome', required=True, help="Input Genome file in FastA format")
    parser.add_argument("-g", "--gff", dest='GFF_File', required=True, help='GFF file from prodigal, option -f "gff"')
    parser.add_argument("-p", "--prot", dest='Protein_File', required=True, help='Protein file from prodigal, option -a')
    parser.add_argument("-o", "--output", dest='Output_File', required=True, help='Output Genbank file')
    parser.add_argument("--annot", dest='Annotation_File', help='Annotation file in tab separated format, please also specify the columns with the gene name and annotation, by default 1 and 3, respectively')
    parser.add_argument("--namecol", dest='Name_Col', help='Column with protein names')
    parser.add_argument("--annotcol", dest='Annotation_Col', help='Column with protein annotations')
    args = parser.parse_args()

    Input_Genome = args.Input_Genome
    GFF_File = args.GFF_File
    Protein_File = args.Protein_File
    Output_File = args.Output_File
    Annotation_File = args.Annotation_File
    Name_Col = args.Name_Col
    Annotation_Col = args.Annotation_Col

    if Name_Col == None:
        Name_Col = 0
    else:
        Name_Col = int(Name_Col) - 1
    if Annotation_Col == None:
        Annotation_Col = 2
    else:
        Annotation_Col = int(Annotation_Col) - 1

    # Run GFF and FastA parser. It creates a temporal GenBank file for the next stepos to read... could be improved.
    print("Processing {}".format(Input_Genome))
    GFF_Fasta_Merger(Input_Genome, GFF_File)

    # Read the tempral GenBank file to include the protein translations
    print("Generating temporal Genbank file...")
    GenBank_Records = SeqIO.read(open("TEMP_GENBANK.gb","r"), "genbank")
    now = datetime.datetime.now()
    #date = '{}-{}-{}'.format(now.day, now.month, now.year)
    date = now.strftime("%d-%b-%Y").upper()
    GenBank_Records.annotations["date"] = str(date)

    # Check if user provided annotations.
    if Annotation_File == None:
        # Pass the record to the function to add the protein translations
        print("Adding protein translations")
        Genbank_File = Protein_2_GenBank(GenBank_Records, Protein_File)
        # Write final Genbank file
        SeqIO.write(Genbank_File, Output_File, "genbank")
        remove("TEMP_GENBANK.gb")
    else:
        # Pass the record to the function to add the protein translations
        print("Adding protein translations")
        Genbank_File = Protein_2_GenBank(GenBank_Records, Protein_File)
        # Pass the new record to the function to add the annotations
        print("Adding protein annotations")
        Genbank_File_Annot = Annotation_2_GenBank(Genbank_File, Annotation_File, Name_Col, Annotation_Col)
        # Write final Genbank file
        SeqIO.write(Genbank_File_Annot, Output_File, "genbank")
        remove("TEMP_GENBANK.gb")


if __name__ == "__main__":
    main()
