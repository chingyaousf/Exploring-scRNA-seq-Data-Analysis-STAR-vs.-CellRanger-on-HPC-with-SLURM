#!/bin/bash
#SBATCH --chdir=../work
#SBATCH --mail-type=ALL
#SBATCH --time=05:00:00
#SBATCH --mem=50000
#SBATCH --nodes=2
#SBATCH --mail-user=chingyao@usf.edu
#SBATCH --job-name=STAR
#SBATCH --output=STAR.out
#SBATCH --partition=rra --qos=rra

## Load necessary modules
module purge
module load apps/fastqc/0.11.5
module load apps/star/2.6.0a
module load apps/samtools/1.3.1

## Create directories
mkdir -p ../FASTQ_QC
mkdir -p ../genomeDir

## Run FastQC
fastqc ../RNA_seq/*.fastq.gz --outdir ../FASTQ_QC

## Define reference and genome directories
genomeFastaFiles=../reference/GCF_000001635.27_GRCm39_genomic.fna
sjdbGTFfile=../reference/genomic.gtf
genomeDir=../genomeDir

## Index the reference genome using STAR
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile

## Define input and output files
genomeDir=../genomeDir
s1red1_R1=../RNA_seq/ERR1211176_1.fastq.gz
s1red1_R2=../RNA_seq/ERR1211176_2.fastq.gz
s1red1_SAM_OUTFILE=../work/Aligned.out.sam
s1red1_BAM_OUTFILE=../work/Aligned.out.bam

## Define the number of processors for sorting
NUM_PROCESSORS=3

## Run STAR alignment and samtools commands
STAR --runThreadN 3 --genomeDir $genomeDir --readFilesCommand gunzip -c --readFilesIn $s1red1_R1 $s1red1_R2 --outSAMtype SAM  
samtools sort -@ $NUM_PROCESSORS -o $s1red1_BAM_OUTFILE $s1red1_SAM_OUTFILE
samtools index $s1red1_BAM_OUTFILE

