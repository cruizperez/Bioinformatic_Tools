#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  2.0
# Date:		 August 12 2019

# Description: This script appends a prefix or suffix to a given list.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""


def Line_Appender(Input_File, Output_File, String, Prefix = False):
    with open(Input_File) as Input, open(Output_File, 'w') as Output:
        if Prefix == True:
            for line in Input:
                line = line.strip()
                Output.write("{}{}\n".format(String,line))
        else:
            for line in Input:
                line = line.strip()
                Output.write("{}{}\n".format(line,String))


################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script appends a suffix or prefix to a list\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input File] -o [Output File] -s [String to add] --prefix [True or False]\n'''
            '''Global mandatory parameters: -i [Input File] -o [Output File] -s [String to add]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='Input_File', action='store', required=True, help='Input file')
    parser.add_argument('-o', '--output', dest='Output_File', action='store', required=True, help='Output file')
    parser.add_argument('-s', '--string', dest='String', action='store', required=True, help='String to add')
    parser.add_argument('--prefix', dest='Prefix', action='store_true', required=False, help='Add as prefix, by default add as suffix')
    args = parser.parse_args()

    Input_File = args.Input_File
    Output_File = args.Output_File
    String = args.String
    Prefix = args.Prefix


    Line_Appender(Input_File, Output_File, String, Prefix)

if __name__ == "__main__":
    main()
