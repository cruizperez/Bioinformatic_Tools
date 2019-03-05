#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import sys, argparse, os, subprocess

"""----------------------------- 1.0 Define Functions -----------------------------"""

def SRA_Downloader(Sra_ID, Output):
    print(Sra_ID, Output)
    Prefix = Sra_ID[0:3]
    print(Prefix)
    Six_chars = Sra_ID[0:6]
    print(Six_chars)
    url = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/" + Prefix + "/" + Six_chars + "/" + Sra_ID + "/" + Sra_ID + ".sra"
    print(url)
    output = "~/scratch/" + Sra_ID + "/" + Sra_ID
    print(output)
    ssh_key = "/nv/hmicro1/cruizperez3/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    print(ssh_key)
    subprocess.call(['ascp', '-i', ssh_key, '-k', '1', '-T', '-l300', url, Output])


def main():
    parser = argparse.ArgumentParser(description='''Given a Blast output and a FastA file, determines which sequences had good matches and retrieves
                                                    either those with good matches or those without for further processing'''
                                    'Global mandatory parameters: [FastA_File] [Blast_File] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-s', '--sra', dest='Sra_ID', action='store', required=True, help='SRA Identifier to download')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output file')
    args = parser.parse_args()

    Sra_ID = args.Sra_ID
    Output_File = args.Output_File

    SRA_Downloader(Sra_ID, Output_File)

if __name__ == "__main__":
    main()
