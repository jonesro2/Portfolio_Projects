#!/bin/bash
#Author: Robin Jones
#RBIF100:Identifying Potential Crispr Genes in a Cohort
mkdir cohort_precrispr
file1=./cohort_topmotifs/*_topmotifs.fasta 
#loop through each topmotif file,save code name,match GG preceeded by any 21 characters
#output gene matches to precrispr file
for file in $file1; do
  name=$(echo $file |awk -F '_' '{print $1'})
  crispr_site=$(grep -E -B1 ".{20}(.GG)" $file |grep -v ^--)
  echo -e "$crispr_site\n" >> ./cohort_precrispr/"$name"_precrispr.fasta
done
