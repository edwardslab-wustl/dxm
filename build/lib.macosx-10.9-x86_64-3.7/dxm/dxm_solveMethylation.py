#!usr/bin/python3

import sys
import getopt
import numpy as np
import math
import dxm
import dxm.DXMfunctions as DXM
import os

def usage():
    print('\'\'\'')
    print('dxm_SolveMethylation.py identifies subclonal methylation patterns in an input file\n')
    print('Usage: '+sys.argv[0]+' -i <input> -d <outDir> -o <outPref>')
    print('\t-i, --input\tThe sampleFileName')
    print('\t-o, --outPref\tPrefix of output files. Default is \'setNamePlease\'')
    print('\nOPTIONS')
    print('\t-c, --maxCoverage\tThe maximum coverage in the data. Important to set for computational time. Default is 400.')
    print('\t-m, --maxSubpop\tThe maximum number of subpopulations to consider. Computation time scales with number of subpopulations. Default is 2')
    print('\t-f, --fracs\tThe input fractions to use. Default is 0 (our precomputation). Set this value as 1 to recompute fractions from a sample. Set at 2 to solve for user supplied fractions')
    print('\t-u, --userSuppliedFracFile\tIf -f (fracs) is set at 2, this is the user supplied fraction file.') 
    print('\t-b, --binFile\tUser supplied bin file. Must be compatible with freq file (see -q)')
    print('\t-q, --freqFile\tUser supplied freq file. Must be compatible with bin file (see -b)')
    print('\'\'\'')

inFile = ''
outPref = 'dxmSetMethylationNamePlease'
maxCov = 400
maxSubpop = 2
fracs = 0
userFracFile = ''
#binFile = sys.path[0]+'/Bins.txt'
#freqFile = sys.path[0]+'/Freqs.txt'
modulePath = os.path.dirname(DXM.__file__)
binFile = modulePath+'/Bins.txt'
freqFile = modulePath+'/Freqs.txt'

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:m:c:f:u:b:q:",["help","input=","outPref=",'maxSubpop','maxCoverage','fracs','userSuppliedFracFile','binFile=','freqFile='])
except getopt.GetoptError:
    usage()
    sys.exit(2)
if len(sys.argv) <= 1:
    usage()
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h","--help"):
        usage()
        sys.exit(2)
    elif opt in ("-i","--input"):
        inFile = arg
    elif opt in ("-o","--outPref"):
        outPref = arg
    elif opt in ("-m","--maxSubpop"):
        try:
            maxSubpop = int(arg)
        except ValueError:
            print('maximum number of subpopulations should be a positive integer.')
            sys.exit(2)
        if maxSubpop < 2:
            print('maximum number of subpopulations should be an integer >= 2.')
            sys.exit(2)
    elif opt in ("-c","--maxCoverage"):
        try:
            maxCov = int(arg)
        except ValueError:
            print('maximum sequencing coverage should be a positive integer.')
            sys.exit(2)
        if maxCov < 4:
            print('maximum sequencing coverage should be a positive integer.')
            sys.exit(2)
    elif opt in ('-f','--fracs'):
        try:
            fracs = int(arg)
        except ValueError:
            print('fracs value should be 1 or 2')
            sys.exit(2)
        if maxCov < 0:
            print('fracs value should be 1 or 2')
            sys.exit(2)
    elif opt in ('-u','--userSuppliedFracFile'):
        userFracFile = arg
    elif opt in ('-b','--binFile'):
        binFile = arg
    elif opt in ('-q','--freqFile'):
        freqFile = arg

#function defined here for debugging only.
def precompute(maxElement):
    myArray = np.zeros((maxElement+1,1))
    myVal = 0.0
    for i in range(1,maxElement+1):
        myVal += math.log((float(i)))
        myArray[i] = myVal
    return myArray

def main():
    refTable = precompute(maxCov)
    validValues = [0,1,2]
    if fracs not in validValues:
        print('fracs value is set incorrectly. If set, should be 1 or 2.')
        sys.exit(2)
    referenceFracs = []
    if fracs == 1:
        referenceFracs = DXM.getBestFractions(maxSubpop,fracs,inFile)
    else:
        referenceFracs = DXM.getBestFractions(maxSubpop,fracs,userFracFile)
    completion = DXM.runDXMnoIP(inFile,outPref,refTable,maxCov,maxSubpop,referenceFracs,binFile,freqFile)

if __name__ == '__main__':
    main()

