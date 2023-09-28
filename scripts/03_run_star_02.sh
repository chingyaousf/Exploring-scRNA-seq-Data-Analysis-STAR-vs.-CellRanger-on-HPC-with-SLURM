#!/bin/bash
#SBATCH --chdir=../work
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00
#SBATCH --mem=50000
#SBATCH --nodes=3
#SBATCH --mail-user=chingyao@usf.edu
#SBATCH --job-name=03_run_star
#SBATCH --output=03_run_star.out
#SBATCH --partition=rra --qos=rra

##close any open program modules and load the modules you need
module purge
module load apps/star/2.6.0a
module load apps/samtools/1.3.1

genomeDir=../genomeDir
s1red1_R1=../RNA_seq/ERR1211176_1.fastq.gz
s1red1_R2=../RNA_seq/ERR1211176_2.fastq.gz
s1red1_SAM_OUTFILE=../work/Aligned.out.sam
s1red1_BAM_OUTFILE=../work/Aligned.out.bam

##Define the number of processors for sorting
NUM_PROCESSORS=3


STAR --runThreadN 3 --genomeDir $genomeDir --readFilesCommand gunzip -c --readFilesIn $s1red1_R1 $s1red1_R2 --outSAMtype SAM  
samtools sort -@ $NUM_PROCESSORS -o $s1red1_BAM_OUTFILE $s1red1_SAM_OUTFILE
samtools index $s1red1_BAM_OUTFILE






