#!/usr/bin/env bash

#===============================================================================
## initialize conda for sh script
# Some Conda package hooks expect this variable to exist during activation.
#===============================================================================

echo "---------------------------------------------------------"
echo "initiallize the conda in bash"
echo "---------------------------------------------------------"
export CONDA_BACKUP_JAVA_HOME="${CONDA_BACKUP_JAVA_HOME:-}"
eval "$(conda shell.bash hook)"

#===============================================================================
### create enviornement for short reads fastqc and fastp and installing
#===============================================================================

echo "---------------------------------------------------------"
echo "setting up enviornment for fastqc and fastp"
echo "---------------------------------------------------------"

conda create -n 01_short_reads_qc -y
conda activate 01_short_reads_qc

##install fastqc and fastp
conda install bioconda::fastqc -y  # for qc evalution
conda install bioconda::fastp -y   # for qc and trimming

echo "---------------------------------------------------------"
echo "setting up enviornment for multiqc"
echo "---------------------------------------------------------"

#===============================================================================
### create enviornment for short reads for multiqc and installling
#===============================================================================

conda create -n 02_short_reads_multiqc -y
conda activate 02_short_reads_multiqc

## install multiqc
conda install bioconda::multiqc -y

echo "---------------------------------------------------------"
echo "multiqc installed"
echo "---------------------------------------------------------"

#===============================================================================
## creating enviornment and installing star for aligning the transcriptomes
#===============================================================================
echo "---------------------------------------------------------"
echo "creating enviornment for star"
echo "---------------------------------------------------------"

conda create -n 03_star -y
conda activate 03_star

echo "---------------------------------------------------------"
echo "installing the star"
echo "---------------------------------------------------------"

conda install bioconda::star -y

echo "---------------------------------------------------------"
echo " successfully install star"
echo "---------------------------------------------------------"

#===============================================================================
#creating enviornment and installing samtools
#===============================================================================

echo "---------------------------------------------------------"
echo "creating enviornment for samtools"
echo "---------------------------------------------------------"

conda create -n 04_samtools -y
conda activate 04_samtools 

echo "---------------------------------------------------------"
echo "installing samtools"
echo "---------------------------------------------------------"

conda install bioconda::samtools -y

echo "---------------------------------------------------------"
echo "samtools install successfully"
echo "---------------------------------------------------------"

#===============================================================================
#creating enviornment and installing psiclass
#===============================================================================

echo "---------------------------------------------------------"
echo "creating env for psiclass and installing it"
echo "---------------------------------------------------------"
conda create -n 05_psiclass -y
conda activate 05_psiclass

echo "---------------------------------------------------------"
conda install bioconda::psiclass -y

echo "---------------------------------------------------------"
echo "successfully install psiclass"
echo "---------------------------------------------------------"