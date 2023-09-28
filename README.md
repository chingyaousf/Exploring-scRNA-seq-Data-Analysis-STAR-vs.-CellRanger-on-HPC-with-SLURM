# Exploring scRNA-seq Data Analysis: STAR vs. CellRanger on HPC with SLURM

## **Overview:**

Single-cell RNA sequencing (scRNA-seq) has opened up new frontiers in genomics, allowing researchers to study gene expression at the individual cell level. When it comes to analyzing scRNA-seq data, two powerful tools, STAR and CellRanger, often come to the forefront. In this article, we will delve into the differences between these two tools and provide a comprehensive guide on how to run them individually on a High-Performance Computing (HPC) cluster using SLURM. We'll also demonstrate how to create a unified SLURM script for both STAR and CellRanger, streamlining the analysis process for researchers.

## **Understanding STAR and CellRanger:**

#### STAR (Spliced Transcripts Alignment to a Reference)

STAR is renowned for its efficiency in RNA-seq read alignment to a reference genome. While it's not designed exclusively for scRNA-seq, it can be a valuable tool in this context due to its speed and accuracy.

#### CellRanger: A Comprehensive Toolkit

CellRanger, developed by 10x Genomics, is tailor-made for scRNA-seq data analysis. It offers an end-to-end solution, simplifying the entire workflow from read alignment to downstream analysis.

## **Running STAR on HPC with SLURM:**

#### Step 1: Check FastQ Quality Control (01_check_fastq_qc)

Before launching STAR, ensure your input FastQ files are of high quality. You can perform this quality check using FastQC or similar tools. Create a SLURM script to submit the FastQC job to the cluster.

#### Step 2: STAR Indexing (02_star_index)

To facilitate efficient alignment, STAR requires a pre-built index of the reference genome. Create a SLURM script to build the STAR index for your chosen reference genome. This step is typically performed only once unless you switch to a different genome.

#### Step 3: Running STAR for Mapping (03_run_star)

Create a SLURM script to execute the STAR alignment command, specifying input FastQ files, the reference genome index, and relevant parameters. After alignment, continue with quantification steps, such as generating gene expression matrices using tools like featureCounts or HTSeq.

## **Running CellRanger on HPC with SLURM:**

#### Step 1: CellRanger Mkfastq (cellranger_mkfastq)

CellRanger Mkfastq is used for demultiplexing raw sequencing data and generating FASTQ files for each sample. Create a SLURM script to submit the CellRanger Mkfastq job, specifying input BCL files and sample information.

#### Step 2: CellRanger Count (cellranger_count)

CellRanger Count is the heart of CellRanger, performing read alignment, barcode processing, and gene expression quantification in one go. Create a SLURM script to execute the CellRanger Count command, specifying input FASTQ files, the reference genome, and sample information. Customize options as needed for specific chemistry, transcriptome reference, and analysis parameters.

## **Streamlining the Process: A Unified SLURM Script:**

To enhance usability and convenience, consider creating a unified SLURM script that combines essential steps from both STAR and CellRanger. This unified script simplifies the analysis process for you and others in the research community. It might include:

-   **Initialization**: Setting environment variables and paths for both STAR and CellRanger.

-   **Read Alignment and Preprocessing**: Incorporating commands for read alignment and preprocessing.

-   **Data Integration and Analysis**: Adding options for data integration and downstream analysis based on research goals.

With this unified script, we provide an efficient and user-friendly solution for scRNA-seq analysis without the need to switch between multiple tools.

## **Differences in FASTQ and Reference Files:**

In the context of single-cell RNA sequencing (scRNA-seq) analysis, understanding the differences in FASTQ and reference files between STAR and CellRanger is crucial.

### STAR Requirements

-   **Sample FASTQ Files**: STAR requires sample FASTQ files containing the sequencing reads.

-   **Reference Genome**: STAR necessitates a reference genome in the form of a genome index, typically created from a reference FASTA file.

-   **Annotation File (GTF)**: An annotation file in GTF (Gene Transfer Format) format is required. This file provides information about gene and exon coordinates, gene names, and other features.

STAR aligns sequencing reads to the reference genome and annotation, enabling the mapping of reads to specific genes and transcripts.

### CellRanger Requirements (developed by 10x Genomics)

-   **Sample FASTQ Files**: Similar to STAR, CellRanger requires sample FASTQ files containing the sequencing reads.

-   **Transcriptome Reference Package**: Instead of a full genome reference, CellRanger uses a pre-built transcriptome reference package. This package includes both a reference FASTA file (usually a transcriptome) and an associated GTF file (or equivalent) containing transcript annotations.

In contrast to STAR, CellRanger does not align reads to the full genome; it aligns them to a pre-built transcriptome. This transcriptome reference package includes a FASTA file with transcript sequences and a GTF file with annotations specific to the experiment or dataset. This approach is suitable for scRNA-seq data analysis because it focuses on quantifying gene and transcript expression rather than aligning reads to the entire genome.

In summary, both STAR and CellRanger perform read alignment and quantification, but CellRanger simplifies the process for single-cell RNA-seq by using a pre-built transcriptome reference that includes transcript annotations, while STAR requires a separate genome reference and annotation file.

## **Conclusion:**

In summary, both STAR and CellRanger offer powerful solutions for scRNA-seq data analysis, each with its unique strengths. STAR excels in read alignment, while CellRanger provides an end-to-end toolkit. Utilizing SLURM on an HPC cluster enhances the scalability and efficiency of your analysis. Additionally, creating a unified SLURM script simplifies the process for researchers, making scRNA-seq analysis accessible and efficient for the entire research community. **Happy analyzing!**

## Reference:

**STAR:**

<https://broadinstitute.github.io/2019_scWorkshop/processing-scrnaseq-data.html#align-the-reads>

[Single Cell RNA-Sequencing of Pluripotent States Unlocks Modular Transcriptional Variation - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S193459091500418X?via%3Dihub)

<https://github.com/alexdobin/STAR>

**CellRanger:**

<https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_fq>

<https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct>

## Blog:

<https://ssidmarine.wordpress.com/2023/09/27/exploring-scrna-seq-data-analysis-star-vs-cellranger-on-hpc-with-slurm/>

## Access data:

**STAR:**

[E-MTAB-2600 \< ArrayExpress \< BioStudies \< EMBL-EBI](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2600?accession=E-MTAB-2600)

<https://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/600/E-MTAB-2600/Files/E-MTAB-2600.idf.txt>

[ENA Browser (ebi.ac.uk)](https://www.ebi.ac.uk/ena/browser/view/ERR1211176?dataType=RUN)

<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.20/>

**CellRanger:**

<https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz>

<https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv>

<https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz>

## **Input & output data available in the data folder**
