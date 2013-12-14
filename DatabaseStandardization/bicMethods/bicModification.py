'''
Created on Dec 13, 2013

@author: mbeoris
'''
#goal of script is to modify BIC database file to simplify
#and standardize
##ONE MAJOR PROBLEM, does not include references nt for insertions and
##deletions, planning on using brca1seq information (BRCA1 U14680) to 
##acquire nt information


from cDNAtoGenomic import get_key_from_value
from cDNAtoGenomic import cDNA_to_genomic
from cDNAtoGenomic import brcaOne


fin = open("BRCA1_u14680.txt", 'rU')
brca1seq =""

for line in fin:
    if '>'not in line:
        brca1seq += line
        brca1seq = brca1seq.replace('\n', '')

#print brca1seq[185-1]
fin.close()



def bicModification(inputfile, outputfile, chr , parentDict, strand=1):
    fin = open(inputfile, 'rU')
    fnew = open(outputfile, 'w')
    fnew.write('chr , start, end, ref, var, clinsig, references')
    fnew.write('\n')
    
    for line in fin:
        line = line.replace(',', '.')
        text_tokens = line.split('\t')
        if "I" in text_tokens[1] or "Exon" in text_tokens[1]:
            fnew.write("")

        else:
            if 'to' in text_tokens[4]:
                nt = str(text_tokens[4])
                ref = nt [0:1]
                var = nt[-1:]
                end = str(text_tokens[2])
                start = end

            elif 'del' in text_tokens[4]:
                nt = str(text_tokens[4])
                nt = nt.replace('del ','')
                var = brca1seq[int(text_tokens[2])-2]
                ref = var + nt
                start = int(text_tokens[2])
                end = start + 1

            elif 'ins' in text_tokens[4]:
                nt = str(text_tokens[4])
                nt = nt.replace('ins ', '')
                ref = brca1seq[int(text_tokens[2])-1]
                var = ref + nt
                end = str(text_tokens[2])
                start = end
            '''    
            results = get_key_from_value(parentDict, int(end))
            end = cDNA_to_genomic(results, strand)
            
            results = get_key_from_value(parentDict, int(start))
            start = cDNA_to_genomic(results, strand)
'''
            
            newLine = str(chr) + ',' + str(start) + ',' + str(end) + ',' + ref + ',' + var +','+ str(text_tokens[13]).upper() +','+ str(text_tokens[26])
            print newLine
            fnew.write(newLine)
            fnew.write('\n')
            


bicModification('brca1_data.txt', 'brca1_data.csv', 17, brcaOne, -1)
#bicModification('testBRCA1_bic.txt', 'testBRCA1_bicCSV.csv', 17, brcaOne, -1)