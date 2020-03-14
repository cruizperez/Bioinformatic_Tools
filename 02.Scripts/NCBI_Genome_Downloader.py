#!/usr/bin/env python

"""
########################################################################
# Author:      Carlos Ruiz, cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# https://github.com/cruizperez/
# Version:    1.0
# Date:      07 Nov 2019

# Description: Downloads multiple genomes from NCBI RefSeq database
using the Assembly Summary file from NCBI. 


########################################################################
"""

################################################################################
"""---1.0 Import Global Modules---"""
import wget 
import time

################################################################################
"""---2.0 Define Functions---"""

def RefSeq_Downloader(Input_File):
    with open(Input_File, 'r') as Input:
        for line in Input:
            line = line.strip().split("\t")
            name = line[0]
            url = line[19]
            folder = url.split("/")[-1]
            new_url = url + "/" + folder + "_genomic.fna.gz"
            print("\nDownloading from {} into {}.fasta.gz\n".format(new_url, name))
            filename = name + ".fasta.gz"
            try:
                filename = wget.download(new_url, out=filename)
            except:
                time.sleep(60)
                print("\nRetrying to download {}".format(new_url))
        

"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''Downloads multiple genomes from NCBI RefSeq database\n'''
            '''using the Assembly Summary file from NCBI.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input_File]\n'''
            '''Global mandatory parameters: -i [Input_File]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='InputFile', action='store', required=True, help='Input list with times as HH:MM:SS, one per line')
    args = parser.parse_args()

    InputFile = args.InputFile
    
    RefSeq_Downloader(InputFile)

if __name__ == "__main__":
    main()
