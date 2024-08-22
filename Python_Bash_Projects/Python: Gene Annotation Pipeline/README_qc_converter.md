#RBIF109 Final Project - QC, file conversion, annotation for paired read sample
#Author: Robin Jones

This purpose of this script is to 1) assess read quality in fastq files 2) index the hg38 reference genome 3) align sample reads to reference genome 4) convert aligned SAM file to BAM file  5)index BAM alignment file and 6) annotate aligned sequences with gene information using a aligned BAM file and hg38knowngene.gtf file.

Execute the following script by copying or typing the following into the command line. Be sure you are in the same directory/folder where you have the script and any required input files saved.
    python3 rjones_final.py

Input (4 files): 
  1. sample1_R1.fastq, sample1_R2.fastq  (sample read files)
  2. hg38.fa (reference genome file)
  3. hg38.knownGene.gtf (gene annotation file)
  Be sure these are in your current working directory along with your script. 

This script will:
  1. Assess the quality of the reads in any .fastq file in your working directory in 2 ways:
    a) using biopython - will output txt file of quality
    b) using FASTQC - will output html report of quality
  2. Index the reference genome hg38 
  3. Align sample reads from fastq files to reference genome and pipe SAM output file to be converted to sorted BAM file. Final output is alignment.sorted.bam
  4. Index the alignment.sorted.bam file for use in IGV
  5. Using the alignment.sorted.bam file as input along with the hg38.knownGene.gtf file, annotate any overlapping sequence regions between the two files by outputting the common sequence regions and their respective annotation information into a new BED file, annot_seqs.bed.    
  

Output locations:
    a. sample1_R1_quality.txt and sample1_R2_quality.txt
    b. sample1_R1_fastqc and sample1_R2_fastqc html files
    c. alignment.sorted.bam and alignment.sorted.bai
    d. annot_seqs
 will all be located in the directory that you ran the script from.  