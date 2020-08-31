#!usr/bin/env python

import numpy as np
cimport numpy as np
DTYPE = np.int
ETYPE = np.double
ctypedef np.double_t ETYPE_t
ctypedef np.int_t DTYPE_t

import math
from timeit import default_timer as timer
from collections import defaultdict

from cython.view cimport array as cvarray
import sys
import random

# helper functions
def load_transitionMatrix(binFile,freqFile):
    bins = []
    freqs = []
    binCount = 0
    freqCount = 0
    with open(binFile,'r') as INPUT:
        for line in INPUT:
            data = line.strip()
            bins.append( float(data) )
            binCount += 1
    with open(freqFile,'r') as INPUT:
        for line in INPUT:
            data = line.strip().split()
            localTup = ( float(data[0]), float(data[1]) )
            freqs.append(localTup)
            freqCount += 1
    if freqCount % 2 == 1:
        print("When loading transition matrix, need even number of frequency lines")
        sys.exit()
    checkVal = 2*(binCount-1)
    if checkVal - freqCount != 0:
        print("Bins and Frequencies for Transition Matrix are not compatible.  2*(numBins-1) should be numFreqs. ")
        sys.exit()
    collapsedFreqs = np.zeros((binCount,4))
    for i in range(0,binCount-1):
        collapsedFreqs[i,0] = freqs[2*i][0]
        collapsedFreqs[i,1] = freqs[2*i][1]
        collapsedFreqs[i,2] = freqs[2*i+1][0]
        collapsedFreqs[i,3] = freqs[2*i+1][1]
    return(np.array(bins),collapsedFreqs)

# Processes an input sample. Sample is a tab delimited file with regionName, position, fractional methylation, sequencing coverage
def load_sample(fileName):
    positions = []
    #orig_positions = []
    orig_positions2 = []
    chromosomes = []
    methVals = []
    coverages = []
    gene = []
    indices = []
    helpLine = defaultdict(int)
    with open(fileName,'r') as INPUT:
        for line in INPUT:
            data = line.strip().split() #split on whitespace
            #data = line.strip().split("\s+")
            #chr1	995363	995364	A1BG	1.0	37
            chrom = data[0] # not used though
            localPos = int(data[1])
            localPos2 = int(data[2]) # not used though
            geneName = data[3] # not used though
            localMeth = float(data[4])
            localCov = float(data[5])
            positions.append(localPos)
            #orig_positions.append(localPos)
            orig_positions2.append(localPos2)
            chromosomes.append(chrom)
            methVals.append(localMeth)
            coverages.append(localCov)
            if geneName not in helpLine.keys():
                gene.append(geneName)
            helpLine[geneName]+=1
    startIndex = 0
    for i in range(0,len(gene)):
        geneKey = gene[i]
        dist = helpLine[geneKey]-1
        endIndex = startIndex + dist
        newTuple = (startIndex,endIndex)
        indices.append(newTuple)
        startIndex = endIndex + 1
    #outTuple = ( gene, np.array(methVals), np.array(positions,dtype = np.int32), np.array(coverages,dtype = np.int32), indices )
    outTuple = ( gene, np.array(methVals), np.array(positions,dtype = np.int32), np.array(coverages,dtype = np.int32), indices, orig_positions2, chromosomes)
    return outTuple

# Processes an input sample. Sample is a tab delimited file with regionName, position, fractional methylation, sequencing coverage
def load_sample_old(fileName):
    positions = []
    methVals = []
    coverages = []
    gene = []
    indices = []
    helpLine = defaultdict(int)
    with open(fileName,'r') as INPUT:
        for line in INPUT:
            data = line.strip().split("\t")
            geneName = data[0] # not used though
            localPos = int(data[1])
            localMeth = float(data[2])
            localCov = float(data[3])
            positions.append(localPos)
            methVals.append(localMeth)
            coverages.append(localCov)
            if geneName not in helpLine.keys():
                gene.append(geneName)
            helpLine[geneName]+=1
    startIndex = 0
    for i in range(0,len(gene)):
        geneKey = gene[i]
        dist = helpLine[geneKey]-1
        endIndex = startIndex + dist
        newTuple = (startIndex,endIndex)
        indices.append(newTuple)
        startIndex = endIndex + 1
    outTuple = ( gene, np.array(methVals), np.array(positions,dtype = np.int32), np.array(coverages,dtype = np.int32), indices )
    return outTuple

#############################
# used to translate user-provided fraction file
def translateUserFracs(fracFile,maxSubpop):
    outVals = []
    with open(fracFile,'r') as INPUT:
        for line in INPUT:
            data = line.strip().split("\t")
            goodFrac = [float(entry) for entry in data]
            totVals = sum(goodFrac)
            normVals = [element/totVals for element in goodFrac]
            outVals.append( np.array(normVals) )
    if len(outVals) != maxSubpop:
        print('User did not supply the same number of fractions to consider as maximum subpopulations')
    return outVals

# finds first element in an array that is greater than the input
def findBin(element,dividers):
    myIndex = 0
    for x in range(0,len(dividers)):
        if element > dividers[x]:
            myIndex += 1
        else:
            break
    return myIndex

# evaluate L1 difference between methylation distribution and expected values from a fraction
def calcFracDifference(methCounts,linearCombo):
    difference = 0
    for idx in range(0,len(methCounts)):
        methVal = idx/100+0.005
        weight = methCounts[idx]
        localDiff = [abs(Z-methVal) for Z in linearCombo]
        minimum = min(localDiff)
        difference += weight*minimum
    return difference

# converts value to binary representation. Indices correspond from large to small (akin to reading binary number left to right)
def toBinary(value,numElement):
    myVals = np.zeros(numElement)
    localValue = value
    for a in range(0,numElement):
        largestPower = 2**(numElement-a-1)
        if localValue >= largestPower:
            myVals[a] += 1
            localValue -= largestPower
    return myVals

# Generates set of linear combinations from domain of fractions
# Includes 0 or 1 as possible methylation values to consider.
def genLinearCombo(fractions,binaryVals):
    numElement = len(fractions)
    fracAsNParray = np.zeros(numElement)
    for i in range(0,numElement):
        fracAsNParray[i] = fractions[i]
    maxElement = 2**numElement
    myLC = []
    for j in range(0,maxElement):
        localBinary = binaryVals[j]
        markerVal = sum( np.multiply(localBinary,fracAsNParray) )
        if markerVal != 0 and markerVal != 1:
            myLC.append(markerVal)
    goodLC = sorted(myLC)
    return goodLC

# generates a random point in a simplex to consider
def genRandomStartSimplex( numElements ):
    tempList = [random.random() for i in range(0,numElements-1)]
    tempList.sort()
    tempList.insert(0,0.0)
    tempList.append(1.0)
    tempFrac = np.array(tempList)
    outFrac = np.diff(tempFrac).tolist()
    outFrac = sorted(outFrac,key=float)
    return outFrac

# Finds the empirical methylation distribution of a sample
# Discretizes methylation for 1% increments
def generateMethDist(methDistFile):
    myMethVals = []
    with open(methDistFile,'r') as INPUT:
        for line in INPUT:
            data = line.strip().split("\t")
            meth = float(data[2])
            myMethVals.append(meth)

    bins = [0.01*i for i in range(0,101)]
    methCounts = [0]*100
    for methVal in myMethVals:
        idx = findBin(methVal,bins)
        if idx >= 100:
            idx = 99
        methCounts[idx] += 1
    totMeth = len(myMethVals)
    newMethCounts = [float(x)/float(totMeth) for x in methCounts]
    return newMethCounts

# Solves for most likely fractions of underlying subpopulations
def solveMostLikelyFrac(methCounts,numSubpop):
    bestSol = []
    bestDiff = 10000000
    numTries = 10000
    maxElement = 2**numSubpop
    binaryVals = []
    for j in range(0,maxElement):
        localBinary = toBinary(j,numSubpop)
        binaryVals.append(localBinary)
    for x in range(0,numTries):
        randomFrac = genRandomStartSimplex(numSubpop)
        localCombo = 0
        try:
            localCombo = genLinearCombo(randomFrac,binaryVals)
        except IndexError:
            print(str(randomFrac))
            print(str(binaryVals))
            print('subpopnum '+str(numSubpop))
        localDiff = calcFracDifference(methCounts,localCombo)
        if localDiff < bestDiff:
            bestDiff = localDiff
            bestSol = randomFrac
    return bestSol

# Provides best fraction to use for Viterbi based on user requirements
def getBestFractions(maxSubpopToConsider,userCompute,methDistFile):
    FRACS = []
    if userCompute == 0: 
        FRACS.append(np.array([1.0]))
        FRACS.append(np.array([0.25,0.75]))
        FRACS.append(np.array([0.08,0.21,0.71]))
        FRACS.append(np.array([0.05,0.12,0.26,0.57]))
    elif userCompute == 2:
        FRACS = translateUserFracs(methDistFile,maxSubpopToConsider) # fracFile if users precompute fractions to consider
    elif userCompute == 1:
        methDist = generateMethDist(methDistFile)
        FRACS.append(np.array([1.0]))
        for numSubpop in range(2,maxSubpopToConsider+1):
            goodOutput = solveMostLikelyFrac(methDist,numSubpop)
            FRACS.append(np.array(goodOutput))
    return FRACS


# functions used for Viterbi
#################################################

# finds index of max element in array
cdef int maxIndDouble(double[:] myArray, int size):
    cdef int max_ind = 0
    cdef int i = 0
    for i in range(size):
        if myArray[i] > myArray[max_ind]:
            max_ind = i
    return max_ind

# finds index of the first element greater than the input in the array
cdef int maxMinInd(double[:] myArray, int size, double dist):
    cdef int i = 0
    cdef int output = size - 1
    for i in range(0,size):
        if myArray[i] >= dist:
            output = i
            return output
    return output

# Stirling approximation for factorials, with log
cdef double logStirling(int variable, np.ndarray[ETYPE_t, ndim=2] refTable, int maxSize):
    cdef double tempVal
    cdef double outVal = 0.0
    if (variable == 0):
        outVal = 0
    elif (variable < maxSize):
        outVal = refTable[variable]
    elif (variable >= maxSize):
        tempVal = float(variable)
        outVal = tempVal*math.log(tempVal)-tempVal+0.5*math.log(2*math.pi*tempVal)
    return outVal

# Calculation of constants for beta prior. Alpha and Beta parameters in Beta prior are set as integers
cdef double logMyBeta(int a, int b, np.ndarray[ETYPE_t,ndim=2] refTable, int maxSize):
    cdef double outVal = 0.0
    outVal = logStirling(a-1,refTable,maxSize)+logStirling(b-1,refTable,maxSize)-logStirling(a+b-1,refTable,maxSize)
    return outVal

# Beta-binomial calculation
cdef double logBetaBinom(int k, int n, int a, int b, np.ndarray[ETYPE_t,ndim=2] refTable, int maxSize):
    cdef double outVal = 0.0
    outVal = logStirling(n,refTable,maxSize)-logStirling(k,refTable,maxSize)-logStirling(n-k,refTable,maxSize)+logMyBeta(k+a,n-k+b,refTable,maxSize)-logMyBeta(a,b,refTable,maxSize)
    return outVal

# Given a number of reads per subpopulation, computes beta-binomial probabilities
# incorporates a little Gaussian noise as well in model.

cdef double translateMeth(int methState, int n, int k, np.ndarray[ETYPE_t,ndim=2] refTable, int maxSize, np.ndarray[ETYPE_t,ndim=2] METHARRAY, np.ndarray[ETYPE_t,ndim=2] UNMETHARRAY):
    cdef double outVal = 0.0

    if n < maxSize:
        if methState == 1:
            outVal = METHARRAY[n,k]
        else:
            outVal = UNMETHARRAY[n,k]
    else:
        print('%s reads is greater than maximum expected coverage' %(n))
        sys.exit()
    return outVal

# Considers all read partitions for beta-binomial computation
cdef double multinomial(int methState, double[:] fractions, int numSubpop, int n, int k, np.ndarray[ETYPE_t,ndim=2] refTable, int maxSize, np.ndarray[ETYPE_t,ndim=2] METHARRAY, np.ndarray[ETYPE_t,ndim=2] UNMETHARRAY, np.ndarray[ETYPE_t,ndim=2] CHOICEARRAY):
    cdef int i,j
    cdef double frac1 = 0.0 # represents unmethylated fractions, 1-frac1 represents methylated
    cdef double tuningFactor = 0.0
    cdef double tuningFactorNKonly = logStirling(k,refTable,maxSize)+logStirling(n-k,refTable,maxSize)-logStirling(n,refTable,maxSize)
    cdef double totProb = 0.0
    cdef double logP = 0.0
    cdef double logQ = 0.0
    cdef double tune1 = 0.0
    cdef double tune2 = 0.0
    for i in range(0,len(fractions)):
        if (methState>>i)&1==0:
            frac1 += fractions[i]
    if numSubpop == 1:
        if frac1 == 0.0:
            # all methylated
            totProb += translateMeth(int(1),n,k,refTable,maxSize,METHARRAY,UNMETHARRAY)
        else:
            # all unmethylated
            totProb += translateMeth(int(0),n,k,refTable,maxSize,METHARRAY,UNMETHARRAY)
    else:
        # guard against precision issues
        if frac1 <= 0.0000001:
            totProb += translateMeth(int(1),n,k,refTable,maxSize,METHARRAY,UNMETHARRAY)
        elif frac1 >= 0.9999999:
            totProb += translateMeth(int(0),n,k,refTable,maxSize,METHARRAY,UNMETHARRAY)
        else:
            logP = math.log(frac1)
            logQ = math.log(1-frac1)
            for i in range(0,n+1):
                for j in range(0,k+1):
                    if j <= i and ((k-j) <= (n-i)):
                        if n < maxSize:
                            tuningFactor = CHOICEARRAY[n,i]
                        else:
                            tuningFactor = logStirling(n,refTable,maxSize)-logStirling(n-i,refTable,maxSize)-logStirling(i,refTable,maxSize)
                        tune1 = tuningFactor + float(i)*logP+float(n-i)*logQ
                        # force call of methylation around 50%
                        totProb += math.exp(tune1)*translateMeth(0,i,j,refTable,maxSize,METHARRAY,UNMETHARRAY)*translateMeth(1,n-i,k-j,refTable,maxSize,METHARRAY,UNMETHARRAY)
    return totProb

# converts to a different base-representation (base-10)
cdef int reconvertBase(int[:] result, int base, int numDigit):
    cdef int output = 0
    cdef int i
    for i in range(0,numDigit):
        output += result[i]*(base**i)
    return output

# converts to different base, e.g. binary representation
cdef void baseConversion(int number, int numSubpop, int base, int[:] result):
    cdef int i = 0
    cdef int factor = 0
    cdef int localNum = number
    for i in range(numSubpop):
        factor = int(math.floor( localNum/pow(base,numSubpop-i-1) ))
        result[numSubpop-i-1] = factor
        localNum -= factor*pow(base,numSubpop-i-1)

# viterbi method.
def viterbi( double[:] x0, double[:,:] Ts, double[:] bins, int[:] pos, int[:] coverages, double[:] fracs, int[:,:] s_star, np.ndarray[ETYPE_t,ndim=2] refTable, np.ndarray[ETYPE_t,ndim=2] METHARRAY, np.ndarray[ETYPE_t,ndim=2] UNMETHARRAY, np.ndarray[ETYPE_t,ndim=2] CHOICEARRAY):
    # static values
    cdef int numSubpop = fracs.shape[0]
    cdef int numBins = bins.shape[0]
    cdef int numCpG = pos.shape[0]   
    cdef int numState = pow(2,numSubpop)

    # variables used within the Viterbi algorithm
    cdef double dist,log_em_prob,log_tr_prob,oldVal,newVal
    cdef int bestPriorStateIndex,prev,prev_bin,j_bin,maxind,tempDepth,methReads
    cdef double frac1 = 0.0
    cdef int i,j,k,methBin = 0

    # memoryviews used during Viterbi computation and for storing results
    V = cvarray(shape=(numCpG,numState), itemsize=sizeof(double),format = "d")
    cdef double [:,:] V_mv = V
    V_mv[:,:] = -sys.maxsize
    for i in range(0,numState):
        V_mv[0,i]= 0.0

    traceback = cvarray(shape=(numCpG,numState), itemsize=sizeof(int), format = "i")
    cdef int [:,:] traceback_mv = traceback
    traceback_mv[:,:] = 0

    binIndicesPython = np.zeros((numSubpop),dtype = np.int32)
    cdef int[:] binIndicesPython_mv = binIndicesPython

    transitionMatrix = np.zeros((4,))
    cdef double[:] TM_mv = transitionMatrix

    # Compute Transition and Emission Probabilities for each CpG (Viterbi Construction)
    for i in range(1,numCpG):
        dist = pos[i] - pos[i-1] # incorporates different transition probabilities depending on distance between CpGs
        if dist <= 0:
            TM_mv[:] = 0.5
            if dist == 0:
                print("CpGs are at same point at %s.  Consider filtering." %(i))
        else:
            maxind = maxMinInd(bins,numBins,dist)
            TM_mv[:] = Ts[maxind-1,:]
        # for each possible methylation state (of underlying subpopulations)
        for j in range(0,numState):
            bestPriorStateIndex = 0
            log_em_prob = 0
            tempDepth = int(coverages[i])
            methReads = int(round(tempDepth*x0[i]))
            log_em_prob = math.log(multinomial(j,fracs,numSubpop,tempDepth,methReads,refTable,refTable.shape[0]-1,METHARRAY,UNMETHARRAY,CHOICEARRAY))
            # consider each previous state, finding the most likely previous state that could lead to each considered state
            for prev in range(0,numState):
                log_tr_prob = 0.0
                numTrans = 0
                for k in range(0,numSubpop):
                    prev_bin = ((prev>>k)&1)
                    j_bin = ((j>>k)&1)
                    log_tr_prob += math.log(TM_mv[2*prev_bin+j_bin])
                oldVal = V_mv[i,j]
                newVal = log_em_prob + log_tr_prob + V_mv[i-1,prev]
                if oldVal < newVal:
                    V_mv[i,j] = newVal
                    bestPriorStateIndex = prev
                elif oldVal == newVal:
                    # tie-breaking in Viterbi
                    if random.randint(0,1) == 0:
                        bestPriorStateIndex = prev
            # update the most likely previous state for the current state
            traceback_mv[i,j] = bestPriorStateIndex

    # variables used for Viterbi traceback.
    cdef int outputIndex,localIndex
    cdef int startIndex = maxIndDouble(V_mv[numCpG-1],numState)    
    cdef double viterbiProb = max(V_mv[numCpG-1])/<double>numCpG

    probPath = [0]*numCpG
    probPath[-1] = viterbiProb*numCpG
    baseConversion(startIndex,numSubpop,2,binIndicesPython_mv) 
    s_star[numCpG-1,:] = binIndicesPython_mv[:]

    for i in range(numCpG-2,-1,-1):
        binIndicesPython_mv[:] = s_star[i+1]
        outputIndex = reconvertBase(binIndicesPython_mv,2,numSubpop)
        localIndex = traceback_mv[i+1,outputIndex]
        baseConversion(localIndex,numSubpop,2,binIndicesPython_mv)
        s_star[i,:] = binIndicesPython_mv[:]
        probPath[i] = V_mv[i,localIndex]
          
    return (viterbiProb,probPath)

###########################################

#builds array for all combinatorial assignments for beta-binomial calculation
def calculateChoice(refTable,maxSize):
    start = timer()
    choiceArray = np.zeros((maxSize,maxSize))
    cdef int i,j
    for i in range(0,maxSize):
        for j in range(0,maxSize):
            if j <= i:
                choiceArray[i,j] = logStirling(i,refTable,maxSize)-logStirling(j,refTable,maxSize)-logStirling(i-j,refTable,maxSize)
    end = timer()
    print('building choice table takes: %s' %(end-start))
    return choiceArray

# function to calculate all beta binomial probabilities
# incorporates predefined beta prior
# returns two arrays with (n,k) where index = value
# methArray = methylated.
# unMethArray = unmethylated.
def calculateAllBetaBinoms(refTable,maxSize):
    start = timer()
    cdef int n,k       
    methArray = np.zeros((maxSize,maxSize))
    unMethArray = np.zeros((maxSize,maxSize))
    cdef double[:,:] methArray_mv = methArray
    cdef double[:,:] unMethArray_mv = unMethArray

    for n in range(1,maxSize):
        for k in range(0,maxSize):
            methArray_mv[n,k] = 0.1*math.exp(logBetaBinom(k,n,200,1,refTable,maxSize))+0.9*math.exp(logBetaBinom(k,n,6,1,refTable,maxSize))
            unMethArray_mv[n,k] = 0.61*math.exp(logBetaBinom(k,n,1,100,refTable,maxSize))+0.39*math.exp(logBetaBinom(k,n,1,3,refTable,maxSize))

    end = timer()
    return (methArray,unMethArray)

########################

#for each region, apply DXM
def runDXMgeneNOIPSubpop(double[:] BINS, double[:,:] FREQS, methVals, pos, coverages, int NUMCPG, np.ndarray[ETYPE_t, ndim=2] refTable, np.ndarray[ETYPE_t,ndim=2] METHARRAY, np.ndarray[ETYPE_t, ndim=2] UNMETHARRAY, np.ndarray[ETYPE_t,ndim=2] CHOICEARRAY,int NUMSUBPOP,LOCALFRAC):

    cdef double[:] METHVALS = methVals
    cdef int[:] POS = pos
    cdef int[:] COVERAGES = coverages
    cdef int i = 0
    cdef double logOdds = 0.0

    (predictedStates, logOdds,probPath) = iterateDXMnoIP(BINS, FREQS, methVals, pos, coverages, NUMCPG, NUMSUBPOP, refTable, LOCALFRAC, METHARRAY, UNMETHARRAY, CHOICEARRAY)

    return( (predictedStates,logOdds,probPath) )

def iterateDXMnoIP(BINS, FREQS, newMeth, newPos, newCoverage, numCpG, numSubpop, refTable, algorithmFrac, METHARRAY, UNMETHARRAY, CHOICEARRAY):
    cdef double[:] algorithmFrac_mv = algorithmFrac
    predictedState = np.zeros((numCpG,numSubpop),dtype = np.int32)
    cdef int[:,:] predictedState_mv = predictedState
    cdef double logOdds = 0
    (logOdds,probPath) = viterbi(newMeth, FREQS, BINS, newPos, newCoverage, algorithmFrac_mv, predictedState_mv, refTable,METHARRAY,UNMETHARRAY,CHOICEARRAY)

    return (predictedState,logOdds,probPath)

# runs DXM solver without integer programming
def runDXMnoIP(sampleName, outPref, refTable,maxCoverage,maxNumSubpop,REFFRACS,BINFILE,FREQFILE):
    (mybins, myfreqs) = load_transitionMatrix(BINFILE,FREQFILE)
    cdef double[:] BINS = mybins # should make memory view of this np array
    cdef double[:,:] FREQS = myfreqs # should make memory view of this np array

    #(GENE, methvals, pos, coverages, INDICES) = load_sample(sampleName)
    (GENE, methvals, pos, coverages, INDICES, pos2, chromosomes) = load_sample(sampleName)
    cdef int i=0
    cdef double[:] METHVALS = methvals
    cdef int[:] POS = pos
    cdef int[:] COVERAGES = coverages
    cdef int NUMCPG = 0
    cdef int maxSIZE = maxCoverage
    cdef int completionStatus = 0
    cdef int numSubpop = 1
    cdef int maxSubpop = maxNumSubpop
    (METHARRAY, UNMETHARRAY) = calculateAllBetaBinoms(refTable,maxSIZE)
    CHOICEARRAY = calculateChoice(refTable,maxSIZE)

    SOLUTIONS = {}
    ALLVIT = defaultdict(list)
    HARDGENES = []
    geneByNumSubpop = defaultdict(list)
    start = timer()
    
    for i in range(0,len(GENE)):
        if i % 1000 == 0:
            timepoint = timer()
            print('finished %s regions in %s' %(i+1,timepoint-start))
        geneName = GENE[i]
        (startIndex, endIndex) = INDICES[i]
        newMeth = np.array(METHVALS[startIndex:endIndex+1])
        newPos = np.array(POS[startIndex:endIndex+1], dtype=np.int32)
        orig_pos2 = pos2[startIndex:endIndex+1]
        chrom = chromosomes[startIndex:endIndex+1]
        newCoverage = np.array(COVERAGES[startIndex:endIndex+1], dtype=np.int32)
        NUMCPG = endIndex-startIndex+1

        numSubpop = 1
        keepIterate = 1
        bestOdds = -10000000
        bestStates = np.zeros((0,0))
        bestProbPath = []
        while keepIterate == 1:
            LOCALFRAC = REFFRACS[numSubpop-1]
            (myStates,logOdds,localProbPath) = runDXMgeneNOIPSubpop(BINS, FREQS, newMeth, newPos, newCoverage, NUMCPG,refTable,METHARRAY,UNMETHARRAY,CHOICEARRAY,numSubpop,LOCALFRAC)            
            if len(myStates) == 0:
                HARDGENES.append(geneName)
                keepIterate = 0
            else:
                if logOdds > bestOdds:
                    bestOdds = logOdds
                    bestStates = myStates
                    bestProbPath = localProbPath
                    if numSubpop == maxSubpop:
                        keepIterate = 0
                        numSubpop -= 1
                else: 
                    keepIterate = 0
                    numSubpop -= 2
                    # cancel natural increment and also go one back, since this current number does not improve fit.
                    # don't need to set anything since prevOdds, prevStates, and prevDMRs should be correct
            ALLVIT[geneName].append(bestOdds)
            numSubpop += 1

        #SOLUTIONS[geneName] = (bestStates,newPos,bestOdds,bestProbPath)
        SOLUTIONS[geneName] = (bestStates,newPos,bestOdds,bestProbPath,orig_pos2,chrom)
        geneByNumSubpop[numSubpop].append(geneName)

    for keyVal in geneByNumSubpop.keys():
        OUTPUT = open('%s_reconstructed_%s_subpops.txt' %(outPref,keyVal),'w')
        targetGenesToWrite = geneByNumSubpop[keyVal]
        for x in range(0,len(targetGenesToWrite)):
            myTarget = targetGenesToWrite[x]
            #(myTracks, myPos,logOdds, myProbPath)= SOLUTIONS[myTarget]
            (myTracks, myPos,logOdds, myProbPath, orig_pos2, chromosomes)= SOLUTIONS[myTarget]
            for y in range(0,len(myTracks)):
                tempTrack = [str(elementToStr) for elementToStr in myTracks[y]]
                outTemp = '\t'.join(tempTrack)
                outPos = myPos[y]
                chrom = chromosomes[y]
                #pos = orig_pos[y]
                pos2 = orig_pos2[y]
                #outLine = '%s\t%s\t%s\n' %(myTarget,outPos,outTemp)
                outLine = '%s\t%s\t%s\t%s\t%s\n' %(chrom, outPos, pos2, myTarget,outTemp)
                OUTPUT.write(outLine)
        OUTPUT.close()
    OUTPUT2 = open('%s_allVitProb.txt' %(outPref),'w')
    for gene in ALLVIT:
        allProbsConsidered = ALLVIT[gene]
        toWrite = ['%s' %(value) for value in allProbsConsidered]
        outLine2 = str(gene)+'\t'+'\t'.join(toWrite)+'\n'
        OUTPUT2.write(outLine2)
    OUTPUT2.close()

    completionStatus = 1
    return completionStatus

def main():
    print("Call functions from this module.  Nothing to run in main class.")

if __name__ == '__main__':
    main()

