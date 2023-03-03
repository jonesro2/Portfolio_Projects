Finding and counting nucleotide motifs in a genome
#Author: Robin Jones

This purpose of this script is to identify which motifs exist in a genome,their total number of recurrences and locations in the genome.

Execute the following script inside the week2 directory by copying or typing the following into the command line:
./motifs1.sh 

This script takes a motif list file(interesting_motifs.txt) and a genome fasta file(r_bifella.fasta) as input. It will:
 1) make a new directory called motifs
 2) read through lines of the motif list one at a time and count the total number of each motif across the r.bifella genome
 3) output/append each motif and its respective count to a file called motif_count.txt (sample output "GGGGG":289)
 4) read through lines of the motif list one at a time and locate each motif within the genes of the r.bifella genome 
 5) ouput/append any r.bifella gene ids that contain the motif along with the corresponding sequence to a file named for that motif(ex: ATAG.fasta). This is done for each of 10 motifs and results in 10 output files.
 
Input files for this script can be found in the week2 directory.
1) interesting_motifs.txt
2) r_bifella.fasta

Output locations:
1) motifs_count.txt will be located in the week2 directory
2) 10 fasta files with gene id and motif locations (ex:"GGGGG.fasta") will be located in the motifs directory.
