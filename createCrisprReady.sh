#!/bin/bash
#Author: Robin Jones
#Identifying Potential Crispr Genes in a Cohort
mkdir cohort_topmotifs
motifs=$(cat motif_list.txt)
file2=./exomesCohort/*.fasta
name=$(basename -a -s .fasta ./exomes/*.fasta) #capture just animal name from *.fasta filename
#print&count motif pattern in each cohort fasta file,no newline,export motif with count to temp. name file for each fasta
for file in $file2; do 
  name=$(basename -a -s .fasta $file)
  for motif in $motifs; do
    echo -n "$motif " >>./cohort_topmotifs/$name.txt
    grep -o "$motif" $file |wc -l >> $name.txt #find motifs in each fasta and count them;output to new file  
  done
#sort motif counts in reverse number order(field2),print 3 highest motif counts to temp. new file
#assign each motif to variable and search for variable in fasta file,output matches to topmotifs file
  sort -nr -k2 $name.txt | head -n 3 |cut -d " " -f1 >./cohort_topmotifs/"$name"_count.txt  
  motif1=$(sed -n '1p' < "$name"_count.txt)
  motif2=$(sed -n '2p' < "$name"_count.txt)
  motif3=$(sed -n '3p' < "$name"_count.txt)
  gene=$(grep -B1 "$motif1\|$motif2\|$motif3" $file | grep -v ^--)
  echo -e "$gene\n" >> ./cohort_topmotifs/"$name"_topmotifs.fasta
  rm "$name"_count.txt
  rm "$name".txt
done 



  




   
