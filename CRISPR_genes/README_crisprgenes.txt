#Identifying Potential Crispr Genes in a Cohort
#Author: Robin Jones

This purpose of this set of scripts is to identify CRISPR gene candidates in organisms with sequenced genomes and that are between 20-30mm in size.

Execute these 5 scripts inside the week4 directory in the following order by copying or typing the following into the command line:
 1)  ./copyExomes.sh 
   This script takes a clinical data file (clinical_data.txt) as input. 
   It will identify sequenced organisms that are between 20 and 30mm in size and copy the fasta files for those organisms to a newly created exomesCohort directory.
   Input files for this script can be found in the week4 directory:
     1) clinical_data.txt
     2) *.fasta files in ./week4/exomes directory --these are the files to be copied
   Output location: codename.fasta files will be copied to exomesCohort directory in week4 directory

2)  ./createCrisprReady.sh
  This script takes a motif_list.txt file and exomesCohort *.fasta files as input.It will:
    a)find and count motifs in each organism's genome and output to temporary file(codename.txt). 
    b)find the 3 highest occurring motifs in each temp file and output to *count.txt temp file.
    c)search for each organism's top 3 motifs in it's genome fasta file and output those genes to *_topmotifs.fasta
  Input files for script are found:
    a) motif_list.txt -- week4 directory
    b) *.fasta files -- exomesCohort directory inside the week4 directory
  Output location: *_topmotifs.fasta files will found in the week4 directory **note:temporary files are removed at the end of script

3) ./identifyCrisprSite.sh
  This script takes *_topmotifs.fasta files as input and will find CRISPR sites in each organism's genome.CRISPR sites are "GG" bases preceeded upstream by any 21 base pairs.
  It will then output genes with CRISPR sites to *_precrispr.fasta files 
  Input files for script are found in the week4 directory.
    a) *_topmotifs.fasta 
  Output location: *_precrispr.fasta files will be found in the week4 directory

4) ./editGenome.sh
  This script takes *_precrispr.fasta files as input and will add an "A" nucleotide before the ".GG" sequence in each gene.
  The edited genes are output to a *_postcrispr.fasta file.
  Input files for script are found in the week4 directory.
    a) *_precrispr.fasta 
  Output location: *_postcrispr.fasta files will be found in the week4 directory

5) ./exomeReport.py
  This script takes the clinical_data.txt and *_postcrispr.fasta files as input.  
  It will produce a report(exome_report.txt) that summarizes the descriptive information about each organism with potential CRISPR sites, gives the total possible number of CRISPR genes across the cohort and lists all the genes with potential sites. 
  Input files for script are found in the week4 directory.
    a) *_postcrispr.fasta 
    b) clinical_data.txt
  Output location: exome_report.txt file will be found in the week4 directory
 
 

