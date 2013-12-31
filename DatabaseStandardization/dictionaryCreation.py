'''
Created on Dec 10, 2013

@author: mbeoris
'''

#v2

#the method below creates a dictionary from a BED file that 
#can be used to translate a cDNA position to a genomic position

#filein:
#will be the BED file that reports the exon ranges

#startSite: 
#if the reported c.DNA is +1 from what is should be startSite=2
#if reported c.DNA is correctly  labeled startSite=1
#if the reported c.DNA is -1 from where it should be startSite=1

#strandOrientation: 
#set this =-1 if reverse strand alignment

def dictionaryCreation(filein, startSite, strandOrientation = 1):
    fin = open (filein, 'rU')
    translateDict = {}
    exonStart = startSite
    
    for line in fin: 
        text_tokens = line.split('\t')
        tmp1 = text_tokens[1]
        tmp2 = text_tokens[2]
        if strandOrientation > 0:
            valueRange = int(tmp2)- int(tmp1)
        else:
            valueRange = int(tmp1)- int(tmp2)
                       
        translateDict[int(tmp1)+1] = range(exonStart, exonStart+valueRange)
        exonStart += valueRange
        
    fin.close()
    return translateDict

    

#note that brca2 reference values are off by one nucleotide, so this function
#works well for this table, Need to implement startsite (0,1,2,etc.) in function
#brca2 example: nt reported 281, actual nt 280 
#this is possibly caused by wrong reference reported on IARC database (HSU gene, not NM)

brca2_dictIARC = dictionaryCreation('BRCA2_Bed.txt', 2)


