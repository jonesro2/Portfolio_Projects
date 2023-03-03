#!/bin/bash
#Author: Robin Jones
#RBIF100:Finding and counting nucleotide motifs
file1=$(cat interesting_motifs.txt) 
file2=r_bifella.fasta 
mkdir motifs
#print&count motif pattern in fasta file,no newline,export to motif_count file
#find,save,remove extra "--" and print gene id number/associated motif sequences,export to motif-name file in new directory 
for line in $file1 
do 
  echo -n "$line": >>motif_count.txt 
  grep -o "$line" $file2 |wc -l >> motif_count.txt 
  motif=$(grep -B1 "$line" $file2 | grep -v ^--)
  echo -e "$motif\n" >> ./motifs/$line.fasta 
done 




