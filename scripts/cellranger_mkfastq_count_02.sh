#!/bin/bash
#SBATCH --chdir=/work/c/chingyao/10X_genomic/yard_02
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --mem=50000
#SBATCH --time=12:00:00
#SBATCH -o run.out
#SBATCH -e run.err
#SBATCH --mail-user=yourmail
#SBATCH --job-name=cellranger_mkfastq_count_job
#SBATCH --output=cellranger_mkfastq_count_job.out
#SBATCH --partition=rra --qos=rra

# Note: Since the tar command does not auto-run after downloading reference transcriptome, 
# if the reference transcriptome is not downloaded in advance, this script should only be run in bash instead of sbatch.
# If the reference transcriptome is downloaded in advance, then you can run sbatch.

# Create the target directories if they don't exist
mkdir -p run_cellranger_mkfastq
mkdir -p run_cellranger_count

# Get base call files (BCL) and store in run_cellranger_mkfastq directory.
wget -P run_cellranger_mkfastq https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
wget -P run_cellranger_mkfastq https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
tar -zxvf run_cellranger_mkfastq/cellranger-tiny-bcl-1.2.0.tar.gz -C run_cellranger_mkfastq/

# Get a reference transcriptome from 10X Genomics support site and store in run_cellranger_count directory.
 wget -P run_cellranger_count https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
fi ; tar -zxvf run_cellranger_count/refdata-gex-GRCh38-2020-A.tar.gz -C run_cellranger_count/

# Activate the cellranger module
module purge
module load apps/cellranger/3.0.2

# Demultiplexing Illumina base call files (BCL) and produce FASTQ files to the run_cellranger_mkfastq directory.
cellranger mkfastq --id=tutorial_walk_through \
  --run=run_cellranger_mkfastq/cellranger-tiny-bcl-1.2.0 \
  --csv=run_cellranger_mkfastq/cellranger-tiny-bcl-simple-1.2.0.csv

# Once you have FASTQ files and a reference transcriptome, you are ready to run cellranger count
# Specify output directories for run_cellranger_count.
cellranger count --id=run_count_test_sample \
   --fastqs=tutorial_walk_through/outs/fastq_path/H35KCBCXY/test_sample \
   --sample=test_sample \
   --transcriptome=run_cellranger_count/refdata-gex-GRCh38-2020-A

