#!/bin/bash

# may need to chmod +x setup.sh
#requires the lyrata.fasta.fa, poly_freq.c and Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.vcf.gz 


############# Setting up Conda Environments ############

# adding conda to path
export PATH=/Users/gracekavanagh/miniconda3/condabin/conda:$PATH
source /Users/gracekavanagh/miniconda3/etc/profile.d/conda.sh




############### Installing dependencies #############

sudo apt install dos2unix
dos2unix lyrata.fasta.fa

## Adding the right channels
conda config --add channels r
conda config --add channels bioconda
conda config --add channels conda-forge

# ## gatk environment setup
# maybe include the -y or the tag that automatically chooses yes
conda create -n gatk_env -y
conda activate gatk_env
conda install -c bioconda gatk4 -y
conda install python -y
conda deactivate

## faststructure
conda create -n fast_structure -y
conda activate fast_structure
conda install faststructure -y
conda install python -y
conda update python -y
conda deactivate

## bcftools (and tabix)
conda create -n bcf_env bcftools tabix -y
conda deactivate

## samtools
conda create -n samtools samtools -y
conda deactivate

## R environment
conda create -n R_env -y
conda activate R_env
conda install r-base -y
conda install bioconductor-biocinstaller -y
		#R packages ; CRAN; Bioconductor; github
		# Cran method ->>  conda install r-dplyr 

conda install r-ggplot2 -y
conda install r-gridExtra -y
conda install r-ggrepel -y
conda install r-vcfR -y
conda install r-adegenet -y
conda install r-StAMPP -y
conda install r-MASS -y
conda install r-adegraphics -y
conda install r-dplyr -y
conda install r-tidyverse -y
conda install r-geosphere -y
conda install r-pegas -y
conda install r-ade4 -y
		# biodoncutcor method ->> conda install bioconductor-ballgown

		# github (devtools) method:vvvvvvvvvvvvv
		# need to isntall devtools first ->> conda install r-devtools
		#activate R console ->> R
		# load devtools ->> library('devtools')
		# devtools::install_github('alyssafrazee/RSkittleBrewer')

		# inside R method? ->> VVV
		# activate R console ->> R
		# use install packages ->> install.packages("name")
conda deactivate

echo "Conda envs set up"

############## Creating indexes and dictionaries ###############

# Index the fasta
conda activate samtools
samtools faidx lyrata.fasta.fa
conda deactivate

# bcftools make index of vcf (if needed)
conda activate bcf_env
bcftools index -t Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.vcf.gz
conda deactivate

# activate conda environment
conda activate gatk_env

# create dictionary of the fasta
gatk CreateSequenceDictionary -R lyrata.fasta.fa

# NB: do not need .tbi of the VCF as this is given

conda deactivate 

############## Compiling C Programme #################

gcc poly_freq.c -o poly_freq -lm
