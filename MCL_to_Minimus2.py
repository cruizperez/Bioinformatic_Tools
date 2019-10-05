#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 Oct 04 2019

# Description: This script parses an MCL clustering output, running
# Minimus2 (AMOS package) to overlap contigs.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def MCL_to_List(MCL_Output, Output_Folder=None, Contig_Folder=None):
    """ Parses a MCL output and returns a dictionary with contigs that belong
        to each cluster and writes a FastA file with unclustered contigs. """
    import pathlib as pl
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    # Folder
    Contig_Folder = pl.Path(Contig_Folder)
    Output_Folder = pl.Path(Output_Folder)
    Output_Folder.mkdir(exist_ok=True, parents=True)
    # Singletons
    Unclustered_Sequences = {}
    Unclustered_Fasta = Output_Folder.joinpath("Unclustered_Contigs.fasta")
    if Unclustered_Fasta.exists():
        Unclustered_Fasta.unlink()
    # List of cluster folders
    Cluster_List = []

    with open(MCL_Output) as MCL:
        Cluster = 1
        Unclustered_Contigs = 0
        for line in MCL:
            Output_Cluster = Output_Folder.joinpath("Cluster_{}".format(Cluster))
            Output_Cluster.mkdir(exist_ok=True, parents=True) # Create cluster folder
            line = line.strip().split()
            if len(line) > 1:
                Out_Fasta = Output_Cluster.joinpath("Cluster_{}.genomes".format(Cluster))
                Cluster_List.append(Out_Fasta) # Open Fasta Output
                with Out_Fasta.open('w') as Outfile:
                    for genome in line:
                        genome = pl.Path(genome).name
                        Genome_file = Contig_Folder.joinpath(genome) # Input fasta folder
                        with Genome_file.open('r') as Outfasta:
                            for title, seq in SimpleFastaParser(Outfasta):
                                Outfile.write(">{}\n{}\n".format(title,seq))
                    Cluster += 1
            else:
                Unclust = Unclustered_Fasta.open('a')
                for genome in line:
                    Unclustered_Contigs += 1
                    genome = pl.Path(genome).name
                    Genome_file = Contig_Folder.joinpath(genome) # Input fasta folder
                    with Genome_file.open('r') as Infasta:
                        for title, seq in SimpleFastaParser(Infasta):
                            Unclust.write(">{}\n{}\n".format(title,seq))
                Unclust.close()

    print("There were {} clusters".format(Cluster))
    print("There were {} unclustered contigs.".format(Unclustered_Contigs))
    return Cluster_List

    # ------------------------

def Minimus2(Cluster_list):
    import pathlib as pl
    import subprocess

    for Cluster in Cluster_list:
        toAmos_File = Cluster.with_suffix('.afg') # Output for script toAmos
        subprocess.Popen(["toAmos", "-s", Cluster, "-o", toAmos_File], stdout=subprocess.PIPE)
    with open("{}_Cluster_Rep.fasta".format(Output_Prefix), 'w') as Output:
        for cluster, genomes in Cluster_Dict.items():
            Sizes = []
            for genome in genomes:
                path = pathlib.Path(genome)
                with open(path, 'r') as Fasta:
                    for title, seq in SimpleFastaParser(Fasta):
                        Sizes.append((title, len(seq)))
            Sizes.sort(key=lambda tup: tup[1], reverse=True)
            print("Finding representative in {} clusters".format(len(Sizes)))
            for values in Sizes:
                genome = values[0][0] + Extension
                Seed_genome = pathlib.Path(ContigDir) / genome
                with open(Seed_genome) as Fasta_Genome:
                    for title, seq in SimpleFastaParser(Fasta_Genome):
                        Output.write(">{}\n{}\n".format(title,seq))

# TODO Try to dereplicate genomes...
                # else:
                #     Query_name = value[0] + ".fa"
                #     Query_genome = pathlib.Path(Contig_Folder) / Query_name
                #     runner = nucmer.Runner(Seed_genome, Query_genome, "Temp_Alignment_{}.coord".format(number), min_id = 95)
                #     runner.run()
                #     number+=1
                #     file_reader = coords_file.reader("Temp_Alignment_{}.coord".format(number))
                #     alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hit


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse
    from sys import argv

    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses a MCL output into clutered and unclustered contigs\n'''
            '''Usage: ''' + argv[0] + ''' -m [MCL Output] -p [Output Prefix] -e [Extension of files]\n'''
            '''Global mandatory parameters: -m [MCL Output] -p [Output Prefix]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    parser.add_argument('-m', '--mclFile', dest='MCL_File', action='store', required=True, help='Input MCL file.')
    parser.add_argument('-o', '--outfolder', dest='Output_Folder', action='store', required=True, help='Prefix of output files.')
    parser.add_argument('-c', '--contigfolder', dest='Contig_Folder', action='store', required=True, help='Prefix of output files.')
    parser.add_argument('-e', '--extension', dest='Extension', action='store', required=True, help='Extension of contig files')
    parser.add_argument('-d', '--directory', dest='ContigDir', action='store', required=True, help='Directory where contigs are located.')
    args = parser.parse_args()

    MCL_File = args.MCL_File
    Output_Folder = args.Output_Folder
    Extension = args.Extension
    ContigDir = args.ContigDir
    Contig_Folder = args.Contig_Folder

    Cluster_List = MCL_to_List(MCL_File, Output_Folder, Contig_Folder)
    Cluster_to_Alignment(Cluster_List)

if __name__ == "__main__":
    main()