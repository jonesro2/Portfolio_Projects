#!bin/python 3

#wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
from Bio import SeqIO
import glob

# Assess quality scores using Biopython and output a quality txt file for each read pair file
fastq_files = glob.glob("*.fastq")
for fastq_file in fastq_files:
    output_filename = fastq_file.split(".")[0] + "_quality.txt"  # Create output filename
    with open(output_filename, "w") as output_file:  # Open file for writing
        for record in SeqIO.parse(fastq_file, "fastq"):  # Parse each record in the FASTQ file
            output_file.write("%s %s\n" % (record.id, record.seq))  # Write record ID and sequence to file
            output_file.write(str(record) + "\n")  # Write full record to file
            output_file.write("Quality Scores: %s\n" % record.letter_annotations["phred_quality"])  # Write quality scores to file


#assess quality scores using fastqc and output html file for each read pair file
import glob
import subprocess
import os

fastq_files=glob.glob("*.fastq")

for fastq_file in fastq_files:
   fastqc_cmd = ['fastqc', fastq_file]
   process = subprocess.Popen(fastqc_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#index the hg38 reference genome using bwa tools
ref="hg38.fa"
index_cmd=['bwa', 'index', ref]
process = subprocess.run(index_cmd, check=True)

read1_fastq="sample1_R1.fastq"
read2_fastq="sample1_R2.fastq"
align_cmd = ['bwa', 'mem', '-t', '4', ref, read1_fastq, read2_fastq]
align_process = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Sort and convert the SAM file to BAM using Samtools
samtools_sort_cmd = ['samtools', 'sort', '-@', '4', '-o', 'alignment.sorted.bam', '-']
samtools_sort_process = subprocess.Popen(samtools_sort_cmd, stdin=align_process.stdout, stdout=subprocess.PIPE)
align_process.stdout.close()
samtools_sort_process.communicate()

#index sorted BAM file (generates .bai)
index_bam_cmd = ['samtools', 'index', 'alignment.sorted.bam']
subprocess.run(index_bam_cmd, check=True)


cmds=["bedtools intersect -abam alignment.sorted.bam -b hg38.knownGene.gtf -wa -wb -bed > annot_seqs.bed"]
# Execute the command using subprocess
process = subprocess.Popen(cmds, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
process.wait()




