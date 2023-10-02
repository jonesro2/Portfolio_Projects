#!bin/python/
from Bio.Seq import Seq
import requests
import sys
#Starting with any gene symbol, this class will find the entrez and ensembl gene IDs, then find the DNA sequence matching the ensembl ID and output it to a FASTA file 
class GeneIdHomologies:
    def __init__(self, gene_symbol):
        self.gene_symbol = gene_symbol
        self.mygene_server = 'https://mygene.info/v3'
        self.ensembl_server= 'http://rest.ensembl.org'
    #find Entrez ID using mygene API
    def GetEntrezId(self):
        endpoint1='/query?q={}&human'.format(self.gene_symbol)
        r = requests.get(self.mygene_server+endpoint1)
        r=r.json()
        entrez_gene_id=r['hits'][0]['_id']
        return entrez_gene_id
    #find Ensembl ID using mygene API
    def GetEnsemblId(self):
        endpoint2 ='/gene/{}?fields=ensembl.gene'.format(self.GetEntrezId())
        r = requests.get(self.mygene_server+endpoint2)
        r=r.json()
        ensembl_gene_id=r['ensembl']['gene']
        return ensembl_gene_id
    # find DNA sequence for ensembl gene ID
    def GetDnaSeq(self):
        header = {"Content-Type": "text/x-fasta"}
        endpoint = "/sequence/id/{}?type=genomic;format=fasta".format(self.GetEnsemblId())
        r = requests.get(self.ensembl_server+endpoint, headers = header)
        with open ('{}.fasta'.format(self.gene_symbol), 'w') as symbol:
            symbol.write(r.text)
    #find gene homologies in other organisms using ensembl gene ID and output alphabetical list
    def GetHomologies(self):
        animal_set=set()
        header = {"Content-Type": "application/json"}   
        endpoint = "/homology/id/{}?format=condensed".format(self.GetEnsemblId())
        r = requests.get(self.ensembl_server+endpoint, headers = header)
        decoded=r.json()
        #get list of homologies then loop through the dictionaries in the list to get the species values and make a set
        homologs=decoded['data'][0]['homologies']
        for organisms in homologs:
            for key, value in organisms.items():
                if key=='species': 
                    animal_set.add(value)
        animal_set.discard('homo_sapiens')
        ordered_animal_set=sorted(animal_set)
        with open ('{}_homology_list.txt'.format(self.gene_symbol.lower()),'w') as homology:
            homology.write("This is a list of organisms that share some homology with the human {} gene. \n".format(self.gene_symbol))
            for animal in ordered_animal_set:
                homology.write("{} \n".format(animal))
gene=GeneIdHomologies('MC1R')
gene.GetDnaSeq()
gene.GetHomologies()

#read in DNA sequence from fasta file and save sequence for later use
def read_fasta(file):
    with open ('{}'.format(file),'r') as fasta:
        lines=[line.strip()for line in fasta]
        global seq
        seq="".join(lines[1:])#get rid of header line
read_fasta('MC1R.fasta')

#initialize variables and set start and stop variables for ORF function
start='ATG'
stop=["TAA","TAG","TGA"]
codon=''
orf_index=[]
orf_reads=[]
stop_index=0
orf_sizes={}
longest_orf=''
#find longest ORF in sequence from read_fasta function
def longest_orf_finder():
    #find each 'ATG' and log it's starting index
    for i in range(len(seq)):
        codon=seq[i:i+3]
        if codon==start:
            orf_index.append(i)
    #for each start codon index, make a list of the sequence reads that follow them
    for index in orf_index:
        orf_reads.append(seq[index:len(seq)])
    #loop through each read to find it's stop(if any) and shorten read to stop codon. update orf_reads with new read.
    for number,read in enumerate(orf_reads):
        for i in range(0,len(read),3):
            if read[i:i+3] in stop:
                stop_index=i+3
                orf_reads[number]=read[0:stop_index] 
                break
   # make list of final read lengths, dictionary with seq as key and length as value, get longest reading frame           
    orf_lengths=[len(string)for string in orf_reads]
    orf_sizes={read:length for read,length in zip(orf_reads,orf_lengths) if len(read)>6}
    global longest_orf
    longest_orf=max(orf_sizes, key=orf_sizes.get)
    return longest_orf
longest_orf_finder()

#get amino acid sequence passing longest_orf from longest_ORF_finder to Seq function and output to previously created fasta file
def translate(gene_symbol):
    my_seq=Seq(longest_orf)
    amino_acid_seq=my_seq.translate()
    with open ('{}.fasta'.format(gene_symbol), 'a') as symbol:
        symbol.write("\n-----------------------------\n")
        symbol.write("\n>{} amino acid sequence".format(gene_symbol.lower()))
        symbol.write("\n{} \n".format(amino_acid_seq))
translate('MC1R')
