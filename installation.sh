#!/usr/bin/env bash
echo "---------------------------------------------------------"
echo "initiallize the conda in bash"
echo "---------------------------------------------------------"
## initialize conda for sh script
# Some Conda package hooks expect this variable to exist during activation.
export CONDA_BACKUP_JAVA_HOME="${CONDA_BACKUP_JAVA_HOME:-}"
eval "$(conda shell.bash hook)"

### create enviornement for short reads fastqc and fastp
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
### create enviornment for short reads for multiqc

conda create -n 02_short_reads_multiqc -y
conda activate 02_short_reads_multiqc

## install multiqc
conda install bioconda::multiqc -y

echo "---------------------------------------------------------"
echo "multiqc installed"
echo "---------------------------------------------------------"

## install star for aligning the transcriptomes
echo "---------------------------------------------------------"
echo "creating enviornment for star"
echo "---------------------------------------------------------"

conda create -n 03_star -y
conda activate 03_star

echo "---------------------------------------------------------"
echo "installing the star"
echo "---------------------------------------------------------"

conda install bioconda::star

echo "---------------------------------------------------------"
echo " successfully install star"
echo "---------------------------------------------------------"