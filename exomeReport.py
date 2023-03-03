#!/bin/python3
#import relevant tools for reading in and looping through files;initialized dict. &lists for later
import csv
import glob
exome_data={}
file_list=[]
final_genes=[]
#open file,select only rows that meet criteria,save relevant codenames as keys, list of other data as values 
#access dict k/v and write to new .txt file incorporated into sentences.
with open("clinical_data.txt", 'r') as clinical_data:
  reader=csv.reader(clinical_data, delimiter='\t')
  next(reader) #skips header row
  for row in reader:
    if row[2]>=str(20) and row[2]<=str(30) and row[4]=="Sequenced":
      stats=[row[0],row[1],row[2],row[3],row[4]]
      exome_data[row[5]]=stats
with open('exome_report.txt', 'w') as exome_doc:
  for key, value in exome_data.items():
    summary=("The {} was discovered by {}, is {}mm in size, and came from a {} environment.".format(key,value[0],value[2],value[3]))   
    exome_doc.write(summary+'\n')
# find files in directory and append to list.
#read each file, save each gene header line to list then make unique gene list without []
for files in glob.glob('*_postcrispr.fasta'):
  file_list.append(files)
for file in file_list:
  with open(file, 'r') as text_file:
    lines=text_file.read().split() #strips \n from strings and makes list
    name_lines=[line[1:]for line in lines[0::2]] #makes list of gene names without >
    for gene in name_lines:
      if gene not in final_genes:
        final_genes.append(gene)
        str_final_genes=",".join(final_genes)
#get count of genes and append count and gene list in a sentence to file created above
number_genes=len(final_genes)
with open('exome_report.txt', 'a') as exome_doc:
  exome_doc.write("\nThere are {} potential crispr genes across the cohort. The genes are: {}".format(number_genes,(str_final_genes)))
    




