#!usr/bin/python3

import sys
#import getopt
import os
from collections import defaultdict
from operator import itemgetter
import numpy as np
import dxm
import argparse

def setup_parser_arguments():
    #global parser
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,description='''
Finds iDMRs from DXM-solved solutions.''')
    parser.add_argument('-v','--vitProbFile',required=True,help="DXM output file with viterbi probabilities.(REQUIRED).")
    parser.add_argument('-m','--methylTraceFile',required=True,help="DXM output file with subpopulation methylation traces.")
    parser.add_argument('-o','--outPref', default="dxm_idmrs", help="Prefix of output files. (Default: dxm_idmrs)")
    parser.add_argument('-n','--numCpG',default=4, type=int, help="Minimum number of CpG expected in an i-DMR. (Default: 4)")
    parser.add_argument('-l','--length', default=50, type=int, help="Minimum length (bp) of an i-DMR. (Default: 50)")
    parser.add_argument('-p','--purity',default=0.99, type=float, help="Minimum purity of i-DMR (fraction of CpG that are differential in an i-DMR). (Default: 0.99)") 

    return parser

parser = setup_parser_arguments()
args = parser.parse_args()

#def usage():
    #print('dxm_calliDMR.py finds iDMRs from DXM-solved solutions.\n')
    #print('Usage: '+sys.argv[0]+' -i <inputFile> -o <outPref>')
    #print('\t-v, --vitProbFile\tDXM output file with viterbi probabilities.')
    #print('\t-m, --methylTraceFile\tDXM output file with subpopulation methylation traces.')
    #print('\t-o, --outPref\tPrefix of output files. Default is \'dxm_dmr\'')
    #print('\nOPTIONS')
    #print('\t-n, --numCpG\tMinimum number of CpG expected in an i-DMR. Default is 4')
    #print('\t-l, --length\tMinimum length (bp) of an i-DMR. Default is 50')
    #print('\t-p, --purity\tMinimum purity of i-DMR (fraction of CpG that are differential in an i-DMR). Default is 0.99')

vitProbFile = args.vitProbFile
methylTraceFile = args.methylTraceFile
outPref = args.outPref
lengthThresh = int(args.length)
numcpgThresh = int(args.numCpG)
purityThresh = float(args.purity)

#try:
    #opts, args = getopt.getopt(sys.argv[1:], "hv:m:o:l:n:p:",["help","vitProbFile=","methylTraceFile=","outPref=",'numCpG=','length=','purity='])
#except getopt.GetoptError:
    #usage()
    #sys.exit(2)
#if len(sys.argv) <= 1:
    #usage()
    #sys.exit(2)
#for opt, arg in opts:
    #if opt in ("-h","--help"):
        #usage()
        #sys.exit(2)
    #elif opt in ("-v","--vitProbFile"):
        #vitProbFile = arg
    #elif opt in ("-m","--methylTraceFile"):
        #methylTraceFile = arg
    #elif opt in ("-o","--outPref"):
        #outPref = arg
    #elif opt in ("-l","--length"):
        #lengthThresh = int(arg)
    #elif opt in ('-n','--numCpG'):
        #numcpgThresh = int(arg)
    #elif opt in ("-p","--purity"):
        #purityThresh = float(arg)

def getGeneToConsider(inFileName):
    putativeDMRgene = {}
    with open(inFileName,'r') as INPUT:
        for line in INPUT:
            data = line.strip().split("\t")
            gene = data[0]
            vitProb1 = float(data[1])
            vitProb2 = float(data[2])
            if vitProb2 - vitProb1 > 0:
                putativeDMRgene[gene] = 1
    return putativeDMRgene

def getDXMdata(sampleFile,putativeDMRgene):
    geneMethTuples = defaultdict(list)
    numSubpop = 0
    with open(sampleFile,'r') as INPUT:
        for line in INPUT:
            data = line.strip().split("\t")
            gene = data[0]
            if numSubpop == 0:
                numSubpop = len(data)-2
            if gene in putativeDMRgene:
                pos = int(float(data[1]))
                methSol = data[2:2+numSubpop]
                goodMethSol = [int(float(x)) for x in methSol]
                geneMethTuples[gene].append( (pos,goodMethSol) )
    return geneMethTuples

def findDiff(methStates):
    firstElement = methStates[0]
    diffVal = all(x==firstElement for x in methStates) 
    return diffVal

def calcDXMDMR(geneMethTuples,outputPrefix):
    # determines DMRs detected by DXM
    DXMdmrs = defaultdict(list)

    for gene in geneMethTuples:
        currDMR = []
        for x in range(0,len(geneMethTuples[gene])):
            element = geneMethTuples[gene][x]
            currPos = element[0]
            methList = element[1]
            areAllSame = findDiff(methList)
            if areAllSame:
                if len(currDMR) != 0:
                    totalCpG = float(len(currDMR)+1)
                    diffHelper = [currDMR[x][1] for x in range(0,len(currDMR))]
                    numDiffCPG = float(sum(diffHelper))
                    totalLength = currDMR[-1][0]-currDMR[0][0]
                    # bias towards merging on the right.
                    if numDiffCPG/totalCpG < purityThresh:
                        if totalCpG < numcpgThresh or totalLength < lengthThresh:
                            currDMR = []
                        else:
                            DXMdmrs[gene].append( (currDMR[0][0],currDMR[-1][0]) )
                            currDMR = []
                    else:
                        currDMR.append( (currPos,0) )
            else:
                currDMR.append( (currPos,1) )

            if x == (len(geneMethTuples[gene]) - 1):
                if len(currDMR) != 0:
                    totalCpG = len(currDMR)
                    diffHelper = [currDMR[x][1] for x in range(0,len(currDMR))]
                    numDiffCPG = sum(diffHelper)
                    totalLength = currDMR[-1][0]-currDMR[0][0]

                    if numDiffCPG/totalCpG >= purityThresh and totalCpG >= numcpgThresh and totalLength >= lengthThresh:
                        DXMdmrs[gene].append( (currDMR[0][0], currDMR[-1][0]) )

    # now process your dmrs and write to output
    OUTPUT = open('%s_DXMdmrs.txt' %(outputPrefix),'w')
    for myKey in DXMdmrs:
        for myGoodDMR in DXMdmrs[myKey]:
            outLine = '%s\t%s\t%s\n' %(myKey,myGoodDMR[0],myGoodDMR[1])
            OUTPUT.write(outLine)
    OUTPUT.close()
    return DXMdmrs

def main():
    filteredGene = getGeneToConsider(vitProbFile)
    DXMdata = getDXMdata(methylTraceFile,filteredGene)
    DXMdmrs = calcDXMDMR(DXMdata,outPref)

if __name__ == '__main__':
    main()

