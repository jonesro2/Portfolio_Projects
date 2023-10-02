#RBIF100 Assignment 3 - Identifying SNPs in D.Gorgon Genome
#Author: Robin Jones

This purpose of this script is to identify the SNPs in D.Gorgon that result in the mold color variation found in patient samples. It will also identify the exact variant associated with each patient sample.

Execute the following script inside the week6 directory by copying or typing the following into the command line:
    python3 pipeline.py -f hawkins_pooled_sequences.fastq
 
 This script takes the harrington_clinical_data.txt, hawkins_pooled_sequences.fastq and dgorgon_reference.fa files as input.  
 
 It will:
    a. make a new fastqs directory and bams directory inside of the week6 directory
    b. trim the low quality bases from end of sequences and corresponding scores from quality scores 
    c. sort the reads from the fastq file by patient barcode, trim off the barcode sequence and output the named fastq files to the fastqs directory 
    d. perform sequence alignment between patient sample reads and reference d.gorgon genomes and output a temporary aligned SAM file
    e. convert SAM files to sorted BAM files,index the BAMs and output all BAMs to bams directory, remove temporary .bam and .sam files
    f. identify and count the bases found at each position along the patient sample genome reads
    g. identify the wt base and SNP by comparing wt d.gorgon genome with patient samples that had 2 different bases at a position
    h. produce a report(report.txt) that gives the SNP base, SNP position, color sample and number of reads for each patient.It also reports the WT base, SNP base, SNP position for each color mold.
  
  Input files for script are found:
    a) hawkins_pooled_sequences.fastq -- week6 directory
    b) harrington_clinical_data.txt -- week6 directory
    c) dgorgon_reference.fa -- week6 directory
  
  Output locations: 
    a) patient-specific FASTQ files are found in the fastqs directory in the week6 directory
    b) sorted and indexed BAM files are found in the bams directory in the week6 directory
    c) report.txt file will be found in the week6 directory