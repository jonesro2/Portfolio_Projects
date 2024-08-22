#RBIF100 Assignment 4 - Identifying MC1R homologs
#Author: Robin Jones

This purpose of this script is to identify MC1R gene IDs, gene and protein sequences and homologs in other organisms.

Note: This script can be customized to any gene symbol by changing the class call parameter from MC1R to any gene symbol. You will need to then alter the file name in the read_fasta function call and gene symbol in the translate function call.

Execute the following script inside the week8 directory by copying or typing the following into the command line:
    python3 homologies.py

  This script only requires a gene symbol for a gene of interest as input (ex: MC1R). 

The script will:
    a. identify the Entrez gene ID for your gene symbol using the mygene API
    b. identify the Ensembl gene ID for your gene symbol using the mygene API
    c. find the DNA sequence using the Ensembl API
    d. find the amino acid sequence using BioPython's myseq and translate functions
    e. identify gene homologies in other organisms using the Ensembl API  
    f. output the gene sequence and amino acid sequence to a FASTA file
    g. output the list of organisms with gene homologies to a gene homology list text file

Output locations: 
    a. gene_name.fasta file will be found in the week8 folder
    b. genename_homology_list.txt file will be found in the week8 folder
