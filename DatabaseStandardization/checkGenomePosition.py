'''
Created on Dec 17, 2013

@author: mbeoris
'''

##compare reported chromosomal position to that within a given chromosome
#for this file i will be accessing chr17.fa that I downloaded from UCSC genome
#browswer "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
#in the future should work on accessing this file from script

from Bio import SeqIO

#record = SeqIO.read("chr17.fa", "fasta")
#region = record.seq[41276045:41276049 ]
#print region

#region2 = record.seq[41251850:41251867 ]
#print region2


record = SeqIO.read("chr9.fa", "fasta")
chr9region = record.seq[0]
print chr9region

fin = open("book1.csv", 'rU')
fout = open('book1ResultsB.txt', 'w')
fout.write("chr9")
fout.write('\n')
for line in fin:
    data = line.split(',')
    #print data[1]
    if 'start' not in data[1]:
        #string1 = "The nucleotide at " 
        #string1 += str(data[1])
        #string1 +=" in the HG19 reference:"
        string2 = record.seq[int(data[1])-1:int(data[2])]
        #fout.write(str(string1))
        fout.write('\t')
        fout.write(str(string2))
        fout.write('\n')

fout.close()    
fin.close()
