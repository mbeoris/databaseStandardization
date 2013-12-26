'''
Created on Dec 13, 2013

@author: mbeoris
'''
#v2.1


#goal of script is to modify BIC database file to simplify
#and standardize
##ONE MAJOR PROBLEM, does not include references nt for insertions and
##deletions, planning on using seqname information (BRCA1 U14680) to 
##acquire nt information


##NEED TO TEST AND MODIFY FOR POSITIVE STRAND
##should probably have input variable of where the gene starts 

##for deletions, start position is the nt before deleted region

from cDNAtoGenomic import get_key_from_value
from cDNAtoGenomic import cDNA_to_genomic
from cDNAtoGenomic import brca1Dict
from dictionaryCreation import dictionaryCreation
from Bio.Seq import Seq


def get_extended_gene_region(filein):
    fin = open(filein, 'rU')
    seqname =""
    for line in fin:
        if '>'not in line:
            line = ''.join(i for i in line if not i.isdigit())
            seqname += line
            seqname = seqname.replace('\n', '').replace('\t', '').replace(" ", '').upper()        
    fin.close()
    return seqname


brca1seq = get_extended_gene_region("BRCA1_extended.txt")
brca2seq = get_extended_gene_region("BRCA2_extended.txt")
brca2Dict = dictionaryCreation('BRCA2_Bed.txt', 2)

def bicModification(inputfile, outputfile, chromosome, parentDict, seqname, codingRegLen, genomicStart, strand=1):
    fin = open(inputfile, 'rU')
    fnew = open(outputfile, 'w')
    fnew.write('chr , start, end, ref, var, clinsig, references, exon/intron, occurrence of mutation')
    fnew.write('\n')

    for line in fin:
        exon = 'exon'
        line = line.replace(',', '.')
        text_tokens = line.split('\t')
        
        if "Exon" in text_tokens[1] or text_tokens[2] == '-':
            continue #skip over first row and any row that does not have an index value
        
        else:
            if "+" in text_tokens[2] or "-" in text_tokens[2]:
                intron= text_tokens[2]
                if not intron[0].isdigit():
                    continue
                for i in range(0, len(intron)):
                    if intron[i] == "+" or intron[i] == "-":
                        operatorPos = i
                cDNAref = int(intron[0:operatorPos])
                results = get_key_from_value(parentDict, cDNAref)
                genomicPos = cDNA_to_genomic(results, strand)
                if intron[operatorPos] == "+":
                    text_tokens[2] = genomicPos + (strand * int(intron[(operatorPos+1):]))
                elif intron[operatorPos] == "-":
                    text_tokens[2] = genomicPos - (strand * int(intron[operatorPos+1:]))
                exon = 'intron'

            else: 
                cDNAref = int(text_tokens[2])
                if cDNAref < codingRegLen:
                    results = get_key_from_value(parentDict, cDNAref)
                    genomicPos = cDNA_to_genomic(results, strand)
                    text_tokens[2] = genomicPos
                else:
                    continue
            
            ntchange = text_tokens[4]
            
            if ('ins' in ntchange and 'del' in ntchange) or ('to' in ntchange and ('del' in ntchange or 'ins' in ntchange)):
                continue
            
            if 'dup' in ntchange:
                ntchange = ntchange.replace('dup', 'ins')
            
            if 'to' in ntchange:
                if len(ntchange) > 6:
                    for i in range(0, len(ntchange)):
                        if ntchange[i:i+2] == 'to':
                            toPos = i
                    for i in range(0, len(ntchange)):
                        if ntchange[i] == "A" or ntchange[i] == "T" or ntchange[i] == "G" or ntchange[i] == "C":
                            lastPos = i
                    ref = ntchange[0:toPos]
                    var = ntchange[(toPos+2):lastPos]
                else:
                    nt = str(ntchange)
                    ref = nt [0:1]
                    var = nt[-1:]
                start = end = int(text_tokens[2])
            
            elif 'del' in ntchange:
                hasNumber = 'false'
                for i in range(0, len(ntchange)):
                    if ntchange[i].isdigit():
                        hasNumber = 'true'
                        lastDig = i            
                nt = str(ntchange)
                nt = nt.replace('del ','').replace(" ", "")
                if strand >0:
                    if hasNumber == "false":
                        var = seqname[((int(text_tokens[2])-1)-genomicStart)]
                        ref = var + nt
                        start = int(text_tokens[2])-1
                        end = start + len(nt)
                    else: 
                        nt = int(nt[0:lastDig-3])
                        var = seqname[((int(text_tokens[2])-1) - genomicStart)]
                        ref = var + seqname[((int(text_tokens[2])- genomicStart)):(int(text_tokens[2])- genomicStart) + nt]
                        start = int(text_tokens[2])-1
                        end = start + int(nt)
                    
                else:
                    nt = str(ntchange)
                    nt = nt.replace('del ','').replace(" ", "")
                    if hasNumber == "false":
                        var = seqname[((int(text_tokens[2])- genomicStart)*strand) + len(nt)]
                        ref = nt + var
                        start = int(text_tokens[2]-len(nt))
                        end = start + len(nt)
                    else:
                        nt = int(nt[0:lastDig-3])
                        var = seqname[((int(text_tokens[2])- genomicStart)*strand) + nt]
                        ref = seqname[((int(text_tokens[2])- genomicStart)*strand):((int(text_tokens[2])- 41277486)*strand) + nt] + var
                        start = int(text_tokens[2]-nt)
                        end = start + int(nt)
                    

                                 

            elif 'ins' in ntchange:
                hasNumber = 'false'
                for i in range(0, len(ntchange)):
                    if ntchange[i].isdigit():
                        hasNumber = 'true'
                #we will not know what is being inserted, so we'll pass over this condition
                if hasNumber == 'true':
                    continue
                nt = str(ntchange)
                nt = nt.replace('ins', '').replace(' ','')
                if strand > 0:
                    ref = seqname[(int(text_tokens[2])- genomicStart)]
                    var = ref + nt
                    start = end = int(text_tokens[2])
                else:
                    ref = seqname[((int(text_tokens[2])- genomicStart)*strand)+1]
                    var = nt + ref
                    start = end = int(text_tokens[2]-1)
                    

            if strand < 0:
                ref = Seq(str(ref))
                ref = ref.reverse_complement()
                var = Seq(str(var))
                var = var.reverse_complement()
                
            newLine = str(chromosome) + ',' + str(start).upper() + ',' + str(end).upper() + ',' + str(ref).replace(" ",'') + ',' + str(var).replace(" ",'') \
             +','+ str(text_tokens[13]).upper() +','+ str(text_tokens[26]) + ',' + exon
            fnew.write(newLine)
            fnew.write('\n')
            
    fnew.close()
            

bicModification('brca2_data.txt', 'brca2_data.csv', 13, brca2Dict, brca2seq, 10987,32889517)
bicModification('brca1_data.txt', 'brca1_data.csv', 17, brca1Dict, brca1seq, 5711,41277486, -1)