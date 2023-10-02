#!bin/python3/
import argparse
import csv
import os
import glob
import re
import pysam

os.mkdir("fastqs")
os.mkdir("bams")

class ParseFastQ(object):
    def __init__(self,filePath,headerSymbols=['@','+']):
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols

    def __iter__(self):
        return self 

    def __next__(self):
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        if nones == 4:
            raise StopIteration
        assert trues == 4,\
            "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
            "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        return tuple(elemList)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", required=True, help="Place fastq inside here")
    args = parser.parse_args()
#parse fastqfile and save each object to a list for multiple iterations
    fastqfile = ParseFastQ(args.fastq)
    fastqobjs=[obj for obj in fastqfile]

#initialize variables for trimming functions
indiv_samples={}
updated_fastqobjs={}
#find and trim the low quality bases from end of sequence and scores from quality score
def trim_ends():
    for obj in fastqobjs:
        bad_grades=["DD","DF","FD","FF"]
        indices=[]
        for grade in bad_grades:
            indices.append(obj[3].find(grade))
        updated_fastqobjs[obj[0]]=[obj[1][0:min(indices)],obj[2],obj[3][0:min(indices)]]
trim_ends()

#read in harrington_clinical_data and save name,barcode to dictionary for multiple iterations
with open ("harrington_clinical_data.txt", 'r') as clinical_data:
    reader=csv.reader(clinical_data, delimiter='\t')
    next(reader) #skips header row
    barcode={row[0]:row[2] for row in reader}
      
#make dictionary of names with genes that match barcode for that name and trim the barcode from each gene
def trim_code():   
    for key,value in barcode.items():
      indiv_samples[key]=[]
      for k,v in updated_fastqobjs.items():
        if v[0][0:5]==value:
          indiv_samples[key]+=[k,v[0][5:],v[1],v[2][5:]]   
#make a new file for each person and write their specific genes to file
    for key,value in indiv_samples.items():
      string_fastq_list="\n".join(value)
      with open ("./fastqs/{}_trimmed.fastq".format(key), 'w') as samples:
        samples.write (string_fastq_list)
trim_code()  

#use os.system to work in command line to generate aligned sam files
def alignment():
    os.system("bwa index dgorgon_reference.fa")
    for fastq_files in glob.glob('./fastqs/*_trimmed.fastq'):
        filename=os.path.basename(fastq_files)
        name=(re.findall("(.+)_",filename)[0])
        os.system ("bwa mem dgorgon_reference.fa ./fastqs/{}_trimmed.fastq > {}.sam".format(name,name))
alignment()

#use os.system to work in command line to convert sam files to sorted bam files
def sam_to_bam():
    for sam_file in glob.glob('*.sam'):
        name=(re.findall("(.+)\.",sam_file)[0])
        os.system('samtools view -bS {}.sam > ./bams/{}.bam'.format(name,name))
        os.system('samtools sort -m 100M -o ./bams/{}.sorted.bam ./bams/{}.bam'.format(name,name))
        os.remove('./bams/{}.bam'.format(name))
        os.remove('{}.sam'.format(name))
sam_to_bam()

#initialize final data dictionaries for report
patient_data={}
mold_color={}

#get dgorgon genome for WT bases
with open ('dgorgon_reference.fa','r') as ref:
    first=ref.readline()
    ref_seq=ref.readline()

#get a dict of names with respective mold color from clinical data
with open('harrington_clinical_data.txt','r') as data:
        reader=csv.reader(data, delimiter='\t')
        next(reader) 
        for row in reader:
            patient_data[row[0]]=[]
            patient_data[row[0]].append(row[1].lower())

            
# identify and locate the SNP for each mold color for each patient      
def pileup():
        #index and loop through each bam, add names and colors from pt_mold dict to final patient_data dict record, convert to samfile alignment
        for bam_file in glob.glob('./bams/*.bam'): 
            pt_name=(os.path.basename(bam_file).split('.')[0])
            os.system('samtools index ./bams/{}.sorted.bam'.format(pt_name))
            samfile = pysam.AlignmentFile("./bams/{}.sorted.bam".format(pt_name), "rb")
            #count up the number of each base in current samfile, add to nt_counts dictionary
            for pileupcolumn in samfile.pileup():
            #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                nt_counts={'A':0,'T':0,'G':0,'C':0}
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                    #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        nt_counts[base]+= 1
                    #only keep the base counts if they are something other than 0 or total reads
                    final_nt_counts={k:v for k,v in nt_counts.items() if v!=0 and v!=pileupcolumn.n}
                #add data from each samfile to final patient_data dictionary for report
                var_pos=' '
                var_snp=' '
                var_freq=' '
                wt_base=' '
                if len(final_nt_counts) !=0:
                    var_pos=pileupcolumn.pos
                    patient_data[pt_name].append(pileupcolumn.n)
                    patient_data[pt_name].append(var_pos)  
                for key,value in final_nt_counts.items():
                    value=round((value/pileupcolumn.n)*100,0) #convert base count to percentage
                    if key != ref_seq[var_pos]: #find which of the two bases in dictionary is the WT/SNP by comparing with dgorgon at same position
                        wt_base=ref_seq[var_pos]
                        patient_data[pt_name].append(wt_base)
                        var_snp=key
                        patient_data[pt_name].append(var_snp)
                        var_freq=value
                        patient_data[pt_name].append(var_freq) 
            samfile.close()  
    
pileup()

#make mold_color dictionary to report out unique color info; report out all patient samples
with open ('report.txt','w') as report:
    report.write("\n Results for Individual Patient Samples\n")
    for key, value in patient_data.items():
        mold_color[value[0]]=[value[2],value[3],value[4]]
        report.write("\n {}'s sample was {} mold and produced {} reads. {}% of the reads at postion {} had the mutation {}.".format(key,value[0],value[1],value[5],value[2],value[4]))   
    report.write("\n---------------------------\n")
    report.write("\n D.Gorgon SNP Evaluation\n")
    for key, value in mold_color.items():
        report.write("\n The {} mold is caused by a SNP in d.gorgon at position {}.The SNP is a change from WT base {} to base {}.".format(key,value[0],value[1],value[2]))
    

       
       
      





        