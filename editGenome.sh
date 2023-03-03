#!/bin/bash
#Author: Robin Jones
#RBIF100:Identifying Potential Crispr Genes in a Cohort
file1=./*_precrispr.fasta
#loop through each precrispr file,save code name,find .GG and add A before that pattern(&)
#output edited genes to postcrispr file
for file in $file1; do
  name=$(echo $file |awk -F '_' '{print $1'})
  crispr_edit=$(sed 's/.GG/A&/' $file) 
  echo -e "$crispr_edit\n" >> "$name"_postcrispr.fasta
done