# dxm

## Installation
Requires python3 (3.5+)

### Required Packages:
     numpy
     cython (>=0.24)
We recommend installing the anaconda (or similar) python distribution, which will include these packages.

### To install:
Download the repository, then:

    pip install .

For user-specific installations, run the install command with the --user flag. 

The "deprecated NumPy API" and "import_array" warnings can be safely ignored.

dxm is also available as a Docker container at: https://hub.docker.com/r/edwardslab/dxm


## Quick Start with Example Data

Data for this example can be found in the example_data folder in the installation folder. The methylation data input file should be in a bed-like format:

<chr> <position1> <position2> <regionName> <fractionalMethylation> <coverage>

### Input file format notes:
regionName - please make unique name for each region tested (e.g. gene name)

position1,position2 - e.g., a genomic coordinate. Please provide as integer.

fractionalMethylation - values should be between 0 (fully unmethylated) and 1 (fully methylated)

coverage - sequencing coverage for that position. Please provide as an integer.

All data should be filtered such that coverage is below the maximum expected sequencing coverage, set as the -m flag in dxm_solveMethylation.


### dxm_estimateFracs
Estimate the fractional prevalence of underlying subpopulations (Optional).

Example: 

    dxm_estimateFracs -i sampleInput.bed -k 3 -o testPrevalence

The output is testPrevalence.txt. Each row is the fractional prevalence of a subpopulation, ordered smallest to largest. Note that this utility is INCOMPATIBLE with dxm_solveMethylation, which has its own fractional prevalence solution call.


### dxm_solveMethylation
Deconvolves processed methylation sequencing data.

Example: 

    dxm_solveMethylation -i sampleInput.bed -o testSample

Notes: 

Set the -m flag to the maximum coverage in the input file.  Increasing -m above the maximum coverage will needlessly increase memory usage.

DXM computational time scales with number of subpopulations, number of CpG, and maximum coverage. DXM solved eRRBS and Methyl-Seq samples for 2 subpopulations at average coverage of 60x in ~2 hours.

The outputs of dxm_solveMethylation are:
- testSample_reconstructed_1_subpops.txt  - regions with 1 major profile
- testSample_reconstructed_2_subpops.txt  - regions with 2 major methylation profiles
- testSample_allVitProb.txt  - list of all relative posterior probabilities

These are tab-delimited files. The format for testSample_reconstructed_2_subpops.txt is:
1. region name
2. position
3. methylation state of minor subpopulation
4. methylation state of major subpopulation


### dxm_callIDMR
Call intrasample differentially methylated regions from solved methylation profiles.

Example: 

    dxm_callIDMR -v testSample_allVitProb.txt -m testSample_reconstructed_2_subpops.txt -o putative

The output of dxm_callIDMR is putative_DXMdmrs.txt. Its format is column delimited
	1. region name
	2. start coordinate
	3. end coordinate

Note: if there are multiple putative iDMRs for the same region, they will have the same corresponding region name.


