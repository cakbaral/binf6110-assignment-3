#!/usr/bin/env bash


# Create the environment for trancriptomics

#Install NCBI SRA Toolkit for obtaining data (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64>
tar -vxzf sratoolkit.tar.gz

#Added the following path to ".bashrc"
# export PATH=$PWD/sratoolkit.3.0.0-mac64/bin:$PATH

#Configure SRA toolkit
vdb-config -i

conda create --name metagenomic
conda activate metagenomic
conda install -c bioconda fastqc
conda install -c bioconda multiqc
conda install -c bioconda fastp
conda install -c bioconda kraken2
conda install -c bioconda bracken
conda install -c bioconda kraken-biom

#----------------------------------

# Download human gut metagenomic samples (from Italy)


#Vegan samples

#Turin - vegan (search result #29, subject ID VOV56)
prefetch SRR8146963 -O /home/cakbarally/binf6110/assignment_3/SRR8146963
fasterq-dump SRR8146963

#Bari - vegan (search result #18, subject ID 01BA)
prefetch SRR8146952 -O /home/cakbarally/binf6110/assignment_3/SRR8146952
fasterq-dump SRR8146952

#Parma - vegan (search result #21, subject ID 26PR)
prefetch SRR8146955 -O /home/cakbarally/binf6110/assignment_3/SRR8146955
fasterq-dump SRR8146955


#Omnivore samples

#Turin - omnivore (search result #41, subject ID VOV77)
prefetch SRR8146975 -O /home/cakbarally/binf6110/assignment_3/SRR8146975
fasterq-dump SRR8146975

#Bari - omnivore (search result #22 - subject ID 26PR)
prefetch SRR8146956 -O /home/cakbarally/binf6110/assignment_3/SRR8146956
fasterq-dump SRR8146956

#Parma - omnivore (search result #62 - subject ID 37PR)
prefetch SRR6367588 -O /home/cakbarally/binf6110/assignment_3/SRR6367588
fasterq-dump SRR6367588



#Download Kraken2 Standard Database (2026-02-26, 16 GB)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20260226.tar.gz
tar -xvzf k2_standard_16_GB_20260226.tar.gz


#----------------------------------

# Quality Control

#Run FASTQC on all files (for raw data)
mkdir -p fastqc_results
fastqc *.fastq -o fastqc_results/ -t 12

#Use MultiQC to obtain a combined report (for raw data)
multiqc fastqc_results/ --outdir multiqc_report_raw

#Trim necessary files (_clean.fastq)
for human in SRR8146963 SRR8146952 SRR8146955 SRR8146975 SRR8146956 SRR6367588
do
echo "Trimming $human..."
fastp \
	-i ${human}_1.fastq \
	-I ${human}_2.fastq \
	-o ${human}_1_clean.fastq \
	-O ${human}_2_clean.fastq \
	-h ${human}_fastp.html \
	-j ${human}_fastp.json
done

#Re-run FASTQ and MultiQC on all files (clean data)
mkdir -p fastqc_clean_results
fastqc *_clean.fastq -o fastqc_clean_results/ -t 12
multiqc fastqc_clean_results/ --outdir multiqc_report_clean

#Move cleaned FASTQ files into separate directory
mkdir -p fastq_clean
mv *_clean.fastq fastq_clean
ls fastq_clean


#----------------------------------

# Taxonomic Classification using Kraken2 and Bracken

#All cleaned files were zipped back for Kraken2.
cd fastq_clean
gzip -v *_clean.fastq

#Output directory
mkdir -p taxonomic_classification


#Kraken2

#NOTE:
#This step was completed on the Narval cluster of DRAC. Please refer to "kraken2.sh" for the exact SLURM script used.


#Bracken
cd taxonomic_classification
for human in SRR8146963 SRR8146952 SRR8146955 SRR8146975 SRR8146956 SRR6367588
do
bracken \
	-d /home/cakbarally/binf6110/assignment_3/kraken_db \
	-i ${human}_kraken2.report \
	-o ${human}.bracken \
	-t 8
done

#Convert to BIOM files

#Separate directory for final BIOM files
mkdir -p BIOM_files

kraken-biom SRR8146963_kraken2_bracken_species.report -o BIOM_files/-o BIOM_files/SRR8146963_table.biom
kraken-biom SRR8146952_kraken2_bracken_species.report -o BIOM_files/-o BIOM_files/SRR8146952_table.biom
kraken-biom SRR8146955_kraken2_bracken_species.report -o BIOM_files/-o BIOM_files/SRR8146955_table.biom
kraken-biom SRR8146975_kraken2_bracken_species.report -o BIOM_files/-o BIOM_files/SRR8146975_table.biom
kraken-biom SRR8146956_kraken2_bracken_species.report -o BIOM_files/-o BIOM_files/SRR8146956_table.biom
kraken-biom SRR6367588_kraken2_bracken_species.report -o BIOM_files/-o BIOM_files/SRR6367588_table.biom

for sample in SRR8146963 SRR8146952 SRR8146955 SRR8146975 SRR8146956 SRR6367588
do
echo "Converting $sample to BIOM..."
kraken-biom ${human}_kraken2_bracken_species.report \
	-o BIOM_files/-o BIOM_files/${human}_table.biom
done


#----------------------------------
# END OF SCRIPT
