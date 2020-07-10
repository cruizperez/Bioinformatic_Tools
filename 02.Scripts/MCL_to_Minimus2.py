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
    from pathlib2 import Path
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    # Folder
    Contig_Folder = Path(Contig_Folder)
    Output_Folder = Path(Output_Folder)
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
                with Out_Fasta.open('w', encoding='utf-8') as Outfile:
                    for genome in line:
                        genome = Path(genome).name
                        Genome_file = Contig_Folder.joinpath(genome) # Input fasta folder
                        with Genome_file.open('r') as Outfasta:
                            for title, seq in SimpleFastaParser(Outfasta):
                                Outfile.write(unicode(">{}\n{}\n".format(title,seq)))
                    Cluster += 1
            else:
                Unclust = Unclustered_Fasta.open('a', encoding='utf-8')
                for genome in line:
                    Unclustered_Contigs += 1
                    genome = Path(genome).name
                    Genome_file = Contig_Folder.joinpath(genome) # Input fasta folder
                    with Genome_file.open('r') as Infasta:
                        for title, seq in SimpleFastaParser(Infasta):
                            Unclust.write(unicode(">{}\n{}\n".format(title,seq)))
                Unclust.close()

    print("There were {} clusters".format(Cluster))
    print("There were {} unclustered contigs.".format(Unclustered_Contigs))

    return Cluster_List

    # ------------------------

def Minimus2(Cluster_list):
    import subprocess

    for Cluster in Cluster_list:
        if Cluster.with_suffix(".out").exists():
            print("{} already processed".format(Cluster))
            continue
        toAmos_File = Cluster.with_suffix('.afg') # Output for script toAmos
        Prefix = Cluster.with_suffix('') # Input prefix for Minimus2
        toAmos = subprocess.Popen(["toAmos", "-s", str(Cluster), "-o", str(toAmos_File)], stdout=subprocess.PIPE)
        toAmos.wait()
        try:
            Minimus = subprocess.check_call(["minimus2", str(Prefix), "-D", "OVERLAP=2000", "-D", "MINID=95"])
        except:
            print "------ WARNING -------"
            print Minimus
        if Minimus == 0:
           Done = Cluster.with_suffix(".out")
           Done.touch(exist_ok=True)

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
    parser.add_argument('-o', '--outfolder', dest='Output_Folder', action='store', required=True, help='Folder to store output.')
    parser.add_argument('-c', '--contigfolder', dest='Contig_Folder', action='store', required=True, help='Folder where individual contigs are located.')
    parser.add_argument('--step', dest="Step", action='store', required=False, default=1, type=int, 
                        help='Step to perform; 1 start from scratch, 2 start from minimus2.')
    args = parser.parse_args()

    MCL_File = args.MCL_File
    Output_Folder = args.Output_Folder
    Contig_Folder = args.Contig_Folder
    Step = args.Step

    if Step == 1:
        Cluster_List = MCL_to_List(MCL_File, Output_Folder, Contig_Folder)
        Minimus2(Cluster_List)
    else:
        from pathlib2 import Path
        Output_Folder = Path(Output_Folder)
        Folders = [f for f in Output_Folder.iterdir() if f.is_dir()]
        Files = [f.joinpath(f.name + '.genomes') for f in Folders]
        Minimus2(Files)

if __name__ == "__main__":
    main()
