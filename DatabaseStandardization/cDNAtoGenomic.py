'''
Created on Dec 5, 2013

@author: mbeoris

goal of module to convert DNA change position from cDNA to genomic position
'''
#
#this is the dictionary for the genomic positions of BRCA1's cds
#key = genomic position of first nucleotide in the exon
#value = range of nucleotides for an exon *note that the exon ranges vary between USCS and genbank


brcaOne = {41277386:    range(1,100+1),     41276132:   range(101,199+1),      41267796:    range(200,253+1), 
           41258550:    range(254,331+1),   41256973:   range(332,420+1),      41256278:    range(421,560+1), 
           41251897:    range(561,666+1),   41249306:   range(667,712+1),      41247939:    range(713,789+1), 
           41246877:    range(790,4215+1),  41243049:   range(4216,4304+1),    41234592:    range(4305,4476+1), 
           41228631:    range(4477,4603+1), 41226538:   range(4604,4794+1),    41223255:    range(4795,5105+1), 
           41219712:    range(5106,5193+1), 41215968:   range(5194,5271+1),    41215390:    range(5272,5312+1), 
           41209152:    range(5313,5396+1), 41203134:   range(5397,5451+1),    41201211:    range(5452,5525+1), 
           41199720:    range(5526,5586+1), 41197819:   range(5587,5711+1)}



def get_key_from_value(my_dict, v):
    for key,value in my_dict.items():
        for i in range(0,len(value)):
            if value[i] == v:
                return key, i
    return False


#this method will convert the cDNA position to genomic position
#input variables are the resultList from get_key_from_value and
#whether or not the gene is on the reverse strand or not (-1 if reverse)
def cDNA_to_genomic(resultList, reverseTF = 1):
    cDNA = -1
    if reverseTF < 0:
        #print "reverse"
        cDNA = resultList[0] - resultList[1]
        return cDNA
    else:
        cDNA = resultList[0] + resultList[1]
    
    return cDNA 
