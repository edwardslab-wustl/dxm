# dxm

DXM performs a local deconvolution of DNA methylation data for one or more regions of the genome.

DXM was first developed by Jerry Fong while he was a member of the Edwards Lab. He is now a resident physician at Barnes Jewish Hospital. DXM is currently maintained by the Edwards Lab.

## Installation
Requires python3 (3.5+)

### Required Python Packages:
     numpy
     cython (>=0.24)
We recommend installing the anaconda (or similar) python distribution, which will include these packages.

### To install:
Download the repository, then:

    pip install .

For user-specific installations, run the install command with the --user flag. 

The "deprecated NumPy API" and "import_array" warnings can be safely ignored.

### Docker Container
dxm is also available as a Docker container at: https://hub.docker.com/r/edwardslab/dxm


## Quick Start with Example Data

Data for this example can be found in the example_data folder in the installation folder. The methylation data input file should be in a bed-like format:

\<chr\> \<position1\> \<position2\> \<regionName\> \<fractionalMethylation\> \<coverage\>

### Input file format notes:

chr - chromosome. Note that this field is not used and can be set to anything.

regionName - please make unique name for each region tested (e.g. gene name)

position1,position2 - e.g., a genomic coordinate. Please provide as integer.

fractionalMethylation - values should be between 0 (fully unmethylated) and 1 (fully methylated)

coverage - sequencing coverage for that position. Please provide as an integer.

All data should be filtered such that coverage is below the maximum expected sequencing coverage, set as the -m flag in dxm_solveMethylation.

DXM generates relative coordinates for internal calculations. As such, it does not explicitly utilize chromosome or position2 data, though these columns are required by DXM to be compatible with BED-like files. DXM computes across all CpGs of a given region, and thus, unique region names should be generated for each region of interest. We recommend adding region names with utilities such as the 'intersect' command from [bedtools](https://bedtools.readthedocs.io/en/latest/).

For example, 
1. Convert your methylation data to a tab-delimited BED-like file with these columns: chromosome position1 position2 methylation_level coverage. See examples below.
2. Download the CGI bed file using the [UCSC Genome Browser Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) to a file named cgi.bed
3. Overlap the data and filter the correct columns: 

    bedtools intersect -wo -a methylation.bed -b cgi.bed | awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5;}' > dxm_in.bed

To convert your methylation data to a tab-delimited BED-like file:

If you are using [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/), start with the .cov produced by bismark_methylation_extractor (chr pos1 pos2 methylation meth_cov unmeth_cov). If your file is called methylation.bismark.cov, then you can convert it using:  

    awk '{cov = $5 + $6; print $1"\t"$2"\t"$3"\t"$4"\t"cov;}' methylation.bismark.cov > methylation.bed

If you are using [bsmap](https://code.google.com/archive/p/bsmap/), use the output from the methratio.py. We recomend also using the -g flag in methratio.py to collapse CpGs across strands.  If your file is called methratio.txt, then you can convert it using:

    awk '{if($4 == "CG") { pos2=$2 + 1; print $1"\t"$2"\t"pos2"\t"$5"\t"$6;}}' methratio.txt > methylation.bed

If you are using [biscuit](https://huishenlab.github.io/biscuit/), the bed output from vcf2bed can be used directly instead of the methylation.bed file in the bed intersect command above.  If you use the mergecg command in the biscuit pipeline (recommended), you must first extract the appropriate columns. If your file after merging is called merge.bed then you can convert it using:

    cut -f 1-5 merge.bed > methylation.bed


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

Set the -c flag to the maximum coverage in the input file.  Increasing -c above the maximum coverage will needlessly increase memory usage.

DXM computational time scales with number of subpopulations, number of CpG, and maximum coverage. DXM solved eRRBS and Methyl-Seq samples for 2 subpopulations at average coverage of 60x in ~2 hours.

The outputs of dxm_solveMethylation are:
- testSample_reconstructed_1_subpops.txt  - regions with 1 major profile
- testSample_reconstructed_2_subpops.txt  - regions with 2 major methylation profiles
- testSample_allVitProb.txt  - list of all relative posterior probabilities

These are tab-delimited files. The format for testSample_reconstructed_1_subpops.txt is:
1. chromosome
2. position
3. position2
4. region name
5. methylation state of major (only) subpopulation

The format for testSample_reconstructed_2_subpops.txt is:
1. chromosome
2. position
3. position2
4. region name
5. methylation state of minor subpopulation
6. methylation state of major subpopulation


### dxm_callIDMR
Call intrasample differentially methylated regions from solved methylation profiles.

Example: 

    dxm_callIDMR -v testSample_allVitProb.txt -m testSample_reconstructed_2_subpops.txt -o putative

The output of dxm_callIDMR is putative_DXMdmrs.txt. Its format is tab-delimited:
	1. chromosome
	2. start coordinate
	3. end coordinate
	4. region name.

Note: if there are multiple putative iDMRs for the same region, they will have the same corresponding region name.


