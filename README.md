






# Title: Investigating allopolyploidy in arenosa and lyrata populations

This project aims to study allopolyploidy by comparing the genetic data from arenosa and lyrata datasets. We used a variety of methods to undergo our analyses:

- Principal Component Analysis (PCA): we used PCA to identify different patterns in the genetic variation of arenosa and lyrata populations to visualise their genetic structures.
- FastStructure Analysis: we used FastStructure to identify the population structure of our data.
- Allele Frequency Spectrum (AFS) Analysis: we conduct AFS analysis on different chromosomes to evaluate the occurance of allopolyploidy.

Through these techniques, we aim to understand the occurance of allopolyploidy in arenosa and lyrata populations.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)



## Installation 


###Before running any of the .R, .py, .C, or .sh files, you need to set up Conda environments and install the necessary dependencies. You can do this by running the 'setup.sh' script. You need to ensure to have the following files in your directory:

# 'lyrata.fasta.fa'
# 'poly_freq.c'
# 'Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.vcf.gz'


Also ensure that the 'setup.sh' file is executable by running:

chmod+x setup.sh

then run the script

./setup.sh




## Usage

To run the project, follow these steps:

### Step 1: Generate genlight PCA

- **Script:** `genlight_PCA.R`
- **Description:** Perform PCA analysis using genlight.
- **Instructions:** Run the `genlight_PCA.R` script to generate the PCA.

### Step 2: Filtering

- **Script:** `filtering_with_gatk_for_filtered_VCF.sh`
- **Description:** Filter the VCF to generate separate and combined VCFs for figures.
- **Instructions:** Run the `filtering_with_gatk_for_filtered_VCF.sh` script to perform filtering.

### Step 3: Toumas PCA

- **Script:** `Toumas_PCA.R`
- **Instructions:** Run the `Toumas_PCA.R` script to generate the PCA.



### Step 4: Splits Tree

- **Tool:** [Splitstree]https://github.com/husonlab/splitstree6?tab=readme-ov-file
- **Description:** Splitstree is a tool used for the analysis of phylogenetic networks. In this project, Splitstree is used to visualize the evolutionary relationships between different populations.
- **Instructions:** Download Splitstree from the link and follow the documentation for analysis.

### Step 5: Trimming and Merging Data

- **Scripts:**
  - `trimming_files.py`
  - `merging_files.py`
- **Description:** Trim and merge data files.
- **Instructions:** Run `trimming_files.py` followed by `merging_files.py` to trim and merge the data.

### Step 6: FastStructure

- **Script:** `faststructure_vcf_generation.py`
- **Description:** Generate VCF for FastStructure analysis.
- **Instructions:** Run the `faststructure_vcf_generation.py` script.


### Step 7: Adgenet

- **Script:** `Adgenet_Ana_Version.R`
- **Description:** Perform Adgenet analysis.
- **Instructions:** Run the `Adgenet_Ana_Version.R` script.

### Step 8: Plot Allele Frequencies

- **Scripts:**
  - `allele_frequency.c`
  - `allele_frequency_plots.R`
- **Description:** Plot allele frequencies.
- **Instructions:** Compile and run `allele_frequency.c`, then run `allele_frequency_plots.R` to plot allele frequencies.

### Step 9: Merging AF Values and Wilcoxon Test

- **Script:** `merging_AF_values_and_Wilcoxon_test.R`
- **Description:** Merge AF values for hybrid population and pure populations and compare them using a Wilcoxon test.
- **Instructions:** Run the `merging_AF_values_and_Wilcoxon_test.R` script.

### Files Needed

- `lyrata.fasta.fa`
- `Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.vcf.gz`
- `pops.txt`: Text file with the names of your populations.



## Credits

### Contributors
- Grace Savanagh
- Vannessa Miller
- Samuel Heysmond

### Supervisor
- Professor Levi Yant
- 


## Data Description

### Input Data

The input data required for this project include:

- `lyrata.fasta.fai`: FASTA file containing reference genome sequences
- 'arenosa.fasta.fai': FASTA file containing reference genome sequences
- `Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.vcf.gz`: Variant Call Format (VCF) file containing genetic variants data.
- Vcf_Name.txt: Text file containing the names of VCF files.
- pops.txt: Tab-delimited file listing individuals to use and their populations
- variables.txt: a text file containing the names of files with allele frequency data




### Output Data

The output data generated from this project include:

- PCA plots: Principal Component Analysis plots visualizing the genetic variation of arenosa and lyrata populations.
- FastStructure results: Population structure analysis results.
- Adgenet analysis results: Genetic clustering analysis results.
- Allele frequency plots: Plots visualizing allele frequencies of different populations.

## File Naming Conventions

- All input and output files should follow the specified naming conventions:
  - `lyrata.fasta.fa`
  - `Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.vcf.gz`

Please ensure that input files are named correctly and are present in the project directory before running the scripts.


### External Libraries and Tools
- [FastStructure](https://github.com/rajanil/fastStructure): Anil Raj, Matthew Stephens, and Jonathan K. Pritchard. fastSTRUCTURE: Variational Inference of Population Structure in Large SNP Data Sets, (Genetics) June 2014 197:573-589 [Genetics, Biorxiv]

- [SplitsTree](https://github.com/husonlab/splitstree6?tab=readme-ov-file):Daniel H. Huson and David Bryant, Application of Phylogenetic Networks in Evolutionary Studies, Molecular Biology and Evolution, 23(2):254-267 (2006) https://doi.org/10.1093/molbev/msj030.

- [Adgenet](https://github.com/thibautjombart/adegenet):Jombart, T., & Ahmed, I. (2011). adegenet: a R package for the multivariate analysis of genetic markers. GitHub repository, https://github.com/thibautjombart/adegenet


