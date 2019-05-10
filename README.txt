dxm

Installation, using python3 (3.5+):

Download dxm-0.1.0.tar.gz from source.

     tar -zxvf dxm-0.1.0.tar.gz
     cd dxm-0.1.0
     python3 setupFull.py install

For user-specific installations, run the install command with the --user flag. 

The "deprecated NumPy API" and "import_array" warnings can be safely ignored.

Required Packages:
     numpy
     cython (>=0.24)
We recommend installing anaconda for a python distribution.

Quick Start with Example Data

Module A) Estimate Fractional Prevalence of underlying subpopulations (Optional)

dxm_estimateFracs -i sampleInput.txt -k 3 -o testPrevalence

The output is testPrevalence.txt. Each row is the fractional prevalence of a subpopulation, ordered smallest to largest. Note that this utility is INCOMPATIBLE with dxm_solveMethylation, which has its own fractional prevalence solution call.

Module B) Deconvolve processed methylation sequencing data

dxm_solveMethylation -i sampleInput.txt -o testSample

The outputs of dxm_solveMethylation are 
	testSample_reconstructed_1_subpops.txt  - regions with 1 major profile
	testSample_reconstructed_2_subpops.txt  - regions with 2 major methylation profiles
	testSample_allVitProb.txt  - list of all relative posterior probabilities

These are column delimited files. The format for testSample_reconstructed_2_subpops.txt is
	1) region name
	2) position
	3) methylation state of minor subpopulation
	4) methylation state of major subpopulation

Module C) Call intrasample differentially methylated regions from solved methylation profiles

dxm_callIDMR -v testSample_allVitProb.txt -m testSample_reconstructed_2_subpops.txt -o putative

The output of dxm_callIDMR is putative_DXMdmrs.txt. Its format is column delimited
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
All data should be filtered such that coverage is below the maximum expected sequencing coverage, set as the -m flag in dxm_solveMethylation.

DXM computational time scales with number of subpopulations, number of CpG, and maximum coverage. DXM solved eRRBS and Methyl-Seq samples for 2 subpopulations at average coverage of 60x in ~2 hours.
