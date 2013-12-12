'''
Created on Dec 10, 2013

@author: mbeoris
'''


#the method below creates a dictionary from a BED file
#that can be used to translate a cDNA position to a genomic
#position, *note change strandOrientation to -1 if reverse strand

def dictionaryCreation(filein, strandOrientation = 1):
    fin = open (filein, 'rU')
    translateDict = {}
    
    
    exonStart = 1
    for line in fin: 
        text_tokens = line.split('\t')
        tmp1 = text_tokens[1]
        tmp2 = text_tokens[2]
        if strandOrientation > 0:
            valueRange = int(tmp2)- int(tmp1)
        else:
            valueRange = int(tmp1)- int(tmp2)
                       
        translateDict[int(tmp1)] = range(exonStart, exonStart+valueRange)
        exonStart += valueRange
        
    keys =  translateDict.keys() 
    print len(keys)
    fin.close()
    return translateDict

    
#dictionaryCreation('test.txt') 

#note that brca2 reference values are off by one nucleotide, so this function
#works well for this table, Need to implement start 0 or start 1 in function
#brca2 example: nt reported 281, actual nt 280
brca2_dict = dictionaryCreation('BRCA2_Bed.txt')


