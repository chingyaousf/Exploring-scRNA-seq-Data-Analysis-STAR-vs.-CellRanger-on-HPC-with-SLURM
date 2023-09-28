#!/bin/bash
#SBATCH --chdir=../work
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00
#SBATCH --mem=10000
#SBATCH --nodes=2
#SBATCH --mail-user=youremail
#SBATCH --job-name=FASTQC
#SBATCH --output=fastqc.out
#SBATCH --partition=rra --qos=rra

##Close any open program modules and load the modules you need
module purge
module load apps/fastqc/0.11.5

##Run the fastqc command. Use a wildcard to specify want all of your fastqs as input. Specify output location
fastqc ../RNA_seq/*.fastq.gz --outdir ../FASTQ_QC

