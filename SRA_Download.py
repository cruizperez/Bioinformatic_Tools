#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.2
# Date:		   25 February 2020

# Description: This script downloades a list of IDs from NCBI SRA achive.
# Requires sra-toolkit and recommended Aspera Connect.
########################################################################
"""

################################################################################

"""---1.0 Import Modules---"""

import sys, argparse, subprocess
import multiprocessing
from shutil import which
from pathlib import Path

################################################################################

"""---2.0 Define Functions---"""

def sra_downloader(sra_id, output, ascp_bin=None, id_rsa=None):
    if ascp_bin != None and id_rsa != None:
        ascp_options = '\'' + str(ascp_bin) + '|' + str(id_rsa) + '\''
        subprocess.call(['prefetch', '-p', '1', '-a', ascp_options, "--ascp-options", '"-k 1 -T -l100m"',
                        "-O", output, sra_id])
    else:
        subprocess.call(['prefetch', '-p', '1', "-O", output, sra_id])

def fastq_dump(sra_id, output, clean=True):
    output_folder = Path(output) / Path(sra_id)
    file_path = (Path(output) / Path(sra_id) / sra_id).with_suffix('.sra')
    subprocess.call(['fastq-dump', '--split-files', '-I', '-O', str(output_folder), str(file_path)])
    if clean == True:
        file_path.unlink()


################################################################################
"""---3.0 Main Function---"""

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                    description='''Download a given SRA dataset and converts it into the corresponding\n'''
                                    '''FastQ files, spliting the file into reads 1 and 2, if possible.\n'''
                                    '''Requires at least sra-toolik, recommended with Aspera Connect\n'''
                                    '''Global mandatory parameters: -l [ID File] or -s [ID List]\n'''
                                    '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-l', '--list_ids', dest='list_ids', action='store', required=False,
                        help='File with SRA id(s) to download (one per line).')
    parser.add_argument('-s', '--sra_ids', dest='Sra_IDs', action='store', nargs='+', required=False,
                        help='List of SRA id(s) to download (separated by spaces).')
    parser.add_argument('-o', '--output', dest='output', action='store', required=False,
                        help='Folder to store files. By default ./')
    parser.add_argument('--aspera_bin', dest='ascp_bin', action='store', required=False,
                        help='''Path to aspera connect, not required but sometimes can be faster.'''
                        '''\nIf used requires aspera\'s private-key file name (--id_rsa).''')
    parser.add_argument('--id_rsa', dest='id_rsa', action='store', required=False,
                        help='Path to aspera\'s private-key file. Usage:\n'''
                        '''--id-rsa [pathtoAspera]/connect/etc/asperaweb_id_dsa.openssh''')
    parser.add_argument('--clean', dest='clean', action='store_false', required=False,
                        help='Remove intermediate .sra files. By default True')
    args = parser.parse_args()

    list_ids = args.list_ids
    Sra_IDs = args.Sra_IDs
    output = args.output
    if output == None:
        output = './'
    ascp_bin = args.ascp_bin
    id_rsa = args.id_rsa

    # Parse input list or file
    list_sra_ids = []
    if list_ids != None and Sra_IDs != None:
        sys.exit("Please provide either a space-separated list of IDs OR a file with one ID per line, not both")
    elif Sra_IDs != None:
        list_sra_ids = Sra_IDs
    elif list_ids != None:
        with open(list_ids, 'r') as list_input:
            for line in list_input:
                line = line.strip()
                list_sra_ids.append(line)
    else:
        sys.exit("No input files provided")

    # Check if both aspera executable and private-key are provided
    if ascp_bin != None and id_rsa == None:
        sys.exit('''I assume you want to use Aspera, for this I need both the path to the executable:\n'''
                '''[pathToAspera]/connect/bin/ascp\n'''
                '''and the private-key file location:\n'''
                '''[pathtoAspera]/connect/etc/asperaweb_id_dsa.openssh''')
    elif ascp_bin == None and id_rsa != None:
        sys.exit('''I assume you want to use Aspera, for this I need both the path to the executable:\n'''
                '''[pathToAspera]/connect/bin/ascp\n'''
                '''and the private-key file location:\n'''
                '''[pathtoAspera]/connect/etc/asperaweb_id_dsa.openssh''')

    # Check if prefetch is available
    if which("prefetch") == None:
        sys.exit('''sra-toolkit binaries not found. Make sure they are in your PATH''')

    # Process ids
    for id in list_sra_ids:
        print("Downloading {}...".format(id))
        sra_downloader(id, output, ascp_bin, id_rsa)
        print("Done")
        print("Processing {}...".format(id))
        fastq_dump(id, output)
        print("Done")

if __name__ == "__main__":
    main()
