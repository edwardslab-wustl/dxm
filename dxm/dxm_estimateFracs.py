#!usr/bin/python3

import sys
#import getopt
import numpy as np
import math
import dxm
import dxm.DXMfunctions as DXM
import os
import argparse

def setup_parser_arguments():
    #global parser
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,description='''
Identifies subclonal methylation patternsi in WGBS data for specific regions.''')
    parser.add_argument('-i','--inFile',required=True,help="Input file name (REQUIRED).  Format should be bed-like: chr pos pos2 regionName methylLevel coverage")
    parser.add_argument('-o','--outPref', default="dxm_fracs", help="Prefix of output files. (Default: dxm_fracs)")
    parser.add_argument('-k','--numSubpop',default=2, type=int, help="The number of subpopulations to consider.Must be >=2. (Default: 2)")

    return parser

#def usage():
#    print('\'\'\'')
#    print('dxm_estimateFracs.py finds the most likely prevalence for a sample given a number of subpopulations\n')
#    print('Usage: '+sys.argv[0]+' -i <input> -k <numSubpop> -o <outPref>')
#    print('\t-i, --input\tThe sampleFileName. Sample is assumed to be in tab-delimited format, with methylation as 3rd column from left.')
#    print('\t-k, --numSubpop\tThe number of subpopulations to consider. Default is 1')
#    print('\t-o, --outPref\tThe output prefix of the file. (default: dxm_fracs)')
#    print('\'\'\'')

parser = setup_parser_arguments()
args = parser.parse_args()

inFile = args.inFile
numSubpop = args.numSubpop
outPref = args.outPref
if numSubpop < 2:
    print('Number of subpopulations should be an integer >= 2.')
    sys.exit(2)

#try:
    #opts, args = getopt.getopt(sys.argv[1:], "hi:k:o:",["help","input=","numSubpop=","outPref="])
##except getopt.GetoptError:
    #usage()
    #sys.exit(2)
#if len(sys.argv) <= 1:
    #usage()
    #sys.exit(2)
#for opt, arg in opts:
    #if opt in ("-h","--help"):
        #usage()
        #sys.exit(2)
    #elif opt in ("-i","--input"):
        #inFile = arg
    #elif opt in ("-k","--numSubpop"):
        #try:
            #numSubpop = int(arg)
        #except ValueError:
            #print('Number of subpopulations should be a positive integer.')
            #sys.exit(2)
        #if numSubpop < 2:
            #print('Number of subpopulations should be an integer >= 2.')
            #sys.exit(2)
    #elif opt in ("-o","--outPref"):
        #outPref = arg

def main():
    methCounts = DXM.generateMethDist(inFile)
    prevalences = DXM.solveMostLikelyFrac(methCounts,numSubpop)
    OUTPUT = open('%s_solvedPrevalences.txt' %(outPref),'w')
    for entry in prevalences:
        OUTPUT.write('%s\n' %(entry))
    OUTPUT.close()

if __name__ == '__main__':
    main()

