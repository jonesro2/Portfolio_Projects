#!/bin/bash 
#Author: Robin Jones
#RBIF100:Identifying Potential Crispr Genes in a Cohort
mkdir exomesCohort
#look through tab-delim file for field matches(size,status), save matched code names field to variable
#loop through names,find matching files in exomes dir, copy to new directory
names=$(awk -F '\t' '{if (($3 >= 20) && ($3 <= 30) && ($5 == "Sequenced")) {print $6}}' clinical_data.txt)
for name in $names
do
  cp ./exomes/$name.fasta ./exomesCohort
done

