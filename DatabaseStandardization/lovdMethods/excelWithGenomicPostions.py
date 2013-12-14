'''
Created on Dec 5, 2013

@author: mbeoris
'''

from cDNAtoGenomic import cDNA_to_genomic
from cDNAtoGenomic import get_key_from_value
from cDNAtoGenomic import brcaOne
from dictionaryCreation import brca2_dict
from Bio.Seq import Seq


def excelWithGenomicPositions(inputFile, outputFile, columnWithcDNAPos, parentDict, chromosome, strand=1):
    fin = open(inputFile, 'rU')
    fout = open(outputFile, 'w')
    startCol = 0
    
    for line in fin:
        text_tokens = line.split(',')
        for i in range(0, len(text_tokens)):
            if columnWithcDNAPos in text_tokens[i]:
                startCol = i
                text_tokens[startCol] = "start"
                text_tokens.insert(startCol+1, "end")
                text_tokens[startCol-1] = "chr"
        
        if startCol !=0 and text_tokens[startCol] != "start":
            cDNApos = int(text_tokens[startCol])
            results = get_key_from_value(parentDict, cDNApos)
            newValue = cDNA_to_genomic(results, strand)
            text_tokens[startCol] = newValue
            text_tokens.insert(startCol+1, newValue)
            text_tokens[startCol-1] = chromosome
            if strand < 0:
                ref = Seq(text_tokens[startCol+2])
                ref = ref.complement()
                ref = str(ref)
                text_tokens[startCol+2] = ref

        text_tokens = str(text_tokens[startCol-1:startCol+4]) + "," + str(text_tokens[startCol+6:startCol+9])
        print text_tokens
        newLine = text_tokens.replace('[','').replace(']', '').replace("'", "")
        newLine = newLine + ' \n'
        fout.write(newLine)  
        
    
    fin.close()
    fout.close()

#excelWithGenomicPositions("test.csv", "testB.csv", " BIC DNA change", brcaOne, 17, -1)
excelWithGenomicPositions("LOVD_BRCA1_12.2.13.csv", "LOVD_BRCA1_12.2.13B.csv", "BIC DNA change", brcaOne, 17, -1)

excelWithGenomicPositions("LOVD_BRCA2_12.10.13.csv", "LOVD_BRCA2_12.10.13B.csv", "BIC DNA change", brca2_dict, 13)