'''
Created on Dec 17, 2013

@author: mbeoris
'''

##compare reported chromosomal position to that within a given chromosome
#for this file i will be accessing chr17.fa that I downloaded from UCSC genome
#browswer "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg9/chromosomes/"
#in the future should work on accessing this file from script

from Bio import SeqIO

#record = SeqIO.read("chr17.fa", "fasta")
#region = record.seq[41276045:41276049 ]
#print region

#region2 = record.seq[41251850:41251867 ]
#print region2


'''record = SeqIO.read("chr13.fa", "fasta")
chr13region = record.seq[0]
print chr13region'''

record = SeqIO.read('chr19.fa', 'fasta')
chr19region = record.seq[50016689-1:50016689+6]
print chr19region
chr19region = record.seq[50029266-1:50029267]
print chr19region
chr19region = record.seq[50017461-1:50017461+6]
print chr19region


'''
fin = open("brca2_dataB.csv", 'rU')
fout = open('brca2results.txt', 'w')
fout.write("chr13")
fout.write('\n')
for line in fin:
    data = line.split(',')
    if 'start' not in data[1]:
        string1 = "The nucleotide at " 
        string1 += str(data[1])
        string1 +=" in the HG19 reference:"
        string2 = record.seq[int(data[1])-1:int(data[2])]
        fout.write(str(string1))
        fout.write('\t')
        fout.write(str(string2))
        fout.write('\n')

fout.close()    
fin.close()
'''