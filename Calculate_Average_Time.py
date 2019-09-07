#!/usr/bin/env python

"""
########################################################################
# Author:   Carlos Ruiz
# Email:    cruizperez3@gatech.edu
# Version:  1.0
# Date:     September 07 2019

# Description: This script calculates the average time from a list
# with the format HH:MM:SS (one per line), it can give the average 
# in hours, minutes or seconds, or as HH:MM:SS.
########################################################################
"""

################################################################################
"""---1.0 Define Functions---"""

def Average_Calculator(InputFile):
    import time
    Hours = []
    Minutes = []
    Seconds = []
    with open(InputFile) as Input_List:
        for line in Input_List:
            line = line.strip().split(":")
            Hours.append(int(line[0]))
            Minutes.append(int(line[1]))
            Seconds.append(int(line[2]))
    Hour_2_Seconds = sum(Hours)*60**2
    Minutes_2_Seconds = sum(Minutes)*60
    Total_Seconds = sum(Seconds) + Hour_2_Seconds + Minutes_2_Seconds
    Seconds_Mean = Total_Seconds / len(Seconds)
    Compound = time.strftime('%H:%M:%S', time.gmtime(Seconds_Mean))

    print("Hours\tMinutes\tSeconds\tCompounds")
    print("{}\t{}\t{}\t{}".format((Seconds_Mean/60**2),(Seconds_Mean/60),Seconds_Mean,Compound))


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script calculates the average time from a list\n'''
            '''with the format HH:MM:SS (one per line), it can give the average\n'''
            '''in hours, minutes or seconds, or as HH:MM:SS.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input List]\n'''
            '''Global mandatory parameters: -i [Input List]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='InputFile', action='store', required=True, help='Input list with times as HH:MM:SS, one per line')
    args = parser.parse_args()

    InputFile = args.InputFile
    
    Average_Calculator(InputFile)

if __name__ == "__main__":
    main()