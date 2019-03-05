#!/usr/bin/env python

"""------------------------- 0.0 Import Modules -----------------------------"""

import sys, argparse, os, subprocess

"""----------------------------- 1.0 Define Functions -----------------------------"""

def SRA_Downloader(Sra_ID, Output, SRA_Bin = None):
    Output_Dir = os.path.dirname(Output)
    if os.path.isdir(Output):
        try:
            print(SRA_ID, " file already found. Attempting to convert to FastQ files")
            if SRA_Bin:
                sra_toolkit_exec = SRA_Bin + "/fastq-dump"
                subprocess.call([sra_toolkit_exec, '--split-files', '-I', '-O', Output_Dir, Output])
            else:
                subprocess.call(['fastq-dump', '--split-files', '-I', '-O', Output_Dir, Output])
        except:
            print("An error occured during processing... Attempting to download file again.")
            pass
    else:
        print("Downloading: ", Sra_ID)
        Prefix = Sra_ID[0:3]
        Six_chars = Sra_ID[0:6]
        url = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/" + Prefix + "/" + Six_chars + "/" + Sra_ID + "/" + Sra_ID + ".sra"
        ssh_key = "/nv/hmicro1/cruizperez3/.aspera/connect/etc/asperaweb_id_dsa.openssh"
        subprocess.call(['ascp', '-i', ssh_key, '-k', '1', '-T', '-l300m', url, Output])
        print(Sra_ID, " downloaded. Now converting to FastQ files (Reads 1 and 2 splitted)")
        if SRA_Bin:
            sra_toolkit_exec = SRA_Bin + "/fastq-dump"
            subprocess.call([sra_toolkit_exec, '--split-files', '-I', '-O', Output_Dir, Output])
        else:
            subprocess.call(['fastq-dump', '--split-files', '-I', '-O', Output_Dir, Output])

def main():
    parser = argparse.ArgumentParser(description='''Download a given SRA dataset and converts it into the corresponding
                                                    FastQ files, as reads 1 and 2 separately.'''
                                    'Global mandatory parameters: [SRA_ID] [Output_File]\n'
                                    'Optional Database Parameters: See ' + sys.argv[0] + ' -h')
    parser.add_argument('-s', '--sra', dest='Sra_ID', action='store', required=True, help='SRA Identifier to download')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output file')
    parser.add_argument('--sra_bin', dest='SRA_Bin', action='store', required=False, help='SRA Toolkit Binaries Folder (required if not on PATH)')
    args = parser.parse_args()

    Sra_ID = args.Sra_ID
    Output_File = args.Output_File
    SRA_Bin = args.SRA_Bin

    if SRA_Bin:
        SRA_Downloader(Sra_ID, Output_File, SRA_Bin)
    else:
        SRA_Downloader(Sra_ID, Output_File)

if __name__ == "__main__":
    main()
