dxm

Installation

Using pip and python3:

download dxm-0.1.0-cp35-cp35m-linux_x86_64.whl
python -m pip install dxm-0.1.0-cp35-cp35m-linux_x86_64.whl

For a user-specific install use the --user flag
python -m pip install --user dxm-0.1.0-cp35-cp35m-linux_x86_64.whl

Required Packages

numpy
cython (>=0.24)
We recommend installing anaconda for a python distribution.

Quick Start with example data

1) Deconvolve processed methylation sequencing data

runDXM -i sampleInput.txt -o testSample

The outputs of runDXM are 
testSample_reconstructed_1_subpops.txt  - regions with 1 major profile
testSample_reconstructed_2_subpops.txt  - regions with 2 major methylation profiles
testSample_allVitProb.txt  - list of all relative posterior probabilities

These are column delimited files. The format for testSample_reconstructed_2_subpops.txt is
	1) region name
	2) position
	3) methylation state of minor subpopulation
	4) methylation state of major subpopulation

2) Call intrasample differentially methylated regions from solved methylation profiles

callIDMR -v testSample_allVitProb.txt -m testSample_reconstructed_2_subpops.txt -o putative

The output of callIDMR is putative_DXMdmrs.txt. Its format is column delimited
	1) region name
	2) start coordinate
	3) end coordinate
Note if there are multiple putative iDMR for the same region, they will have the same corresponding region name.

User Notes

The format for the input methylation sequencing data is column delimited.
	1)  region name - please make unique name for each region tested (e.g. gene name)
	2)  position - e.g., a genomic coordinate. Please provide as integer.
	3)  fractional methylation - between 0 (fully unmethylated) and 1 (fully methylated)
	4)  coverage - sequencing coverage for the CpG. Please provide as integer.
All data should be filtered such that coverage is below the maximum expected sequencing coverage, set as the -m flag in runDXM.

DXM uses TBD computational resources.
