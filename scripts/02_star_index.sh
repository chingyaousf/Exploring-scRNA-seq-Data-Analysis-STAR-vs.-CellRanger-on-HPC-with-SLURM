#!/bin/bash
#SBATCH --chdir=../work
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00
#SBATCH --mem=50000
#SBATCH --nodes=2
#SBATCH --mail-user=chingyao@usf.edu
#SBATCH --job-name=02_star_index
#SBATCH --output=02_star_index.out
#SBATCH --partition=rra --qos=rra

##close any open program modules and load the modules you need
module purge
module load apps/star/2.6.0a

genomeFastaFiles=../reference/GCF_000001635.27_GRCm39_genomic.fna

sjdbGTFfile=../reference/genomic.gtf

genomeDir=../genomeDir

##The STAR will index your reference genome allowing faster alignment
##Input is the reference genome and gtf filles 
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile

