'''
Created on Dec 12, 2013

@author: mbeoris
'''
#goal of this module is to find entries from clinvar_YYYYMMDD.vcf file
#input values will be clinvar filename, the newfile will be a condensed version of
#the original clinvar file


def VCF_Modification(inputfile,newfile):
    clinsig =''
    clnbn =''
    clnorigin =''
    refinfo =''
    brca1count =0
    brca2count =0
    fin = open(inputfile, 'rU')
    fout = open(newfile, 'w')

    for line in fin:
        if 'GENEINFO=BRCA1' in line:
            brca1count += 1
        if 'GENEINFO=BRCA2' in line:
            brca2count +=1
        if '#CHROM' in line:
            fout.write("chr" + '\t' +  "start" + '\t' + "end"+ '\t' + "ref" + '\t' +"var"+ '\t'+ "CLINSIG" + '\t'+"CLNBN" + '\t'+ "CLNORIGIN"+ '\t' +"INFO")
            fout.write('\n')
            
        elif "##" not in line:
            line = line.replace(",", " or ")
            text_tokens = line.split('\t')
            
            
            findInfo = (text_tokens[7]).split(';')
            for i in range(0,len(findInfo)):
                if "CLNSIG" in findInfo[i]:
                    clinsig = findInfo[i]
                    clinsig = clinsig.replace("CLNSIG=", '')
                elif "CLNDBN" in findInfo[i]:
                    clnbn =  findInfo[i]
                    clnbn = clnbn.replace('CLNDBN=', '')
                elif "CLNORIGIN" in findInfo[i]:
                    clnorigin = findInfo[i]
                    clnorigin = clnorigin.replace("CLNORIGIN=",'')
                    if clnorigin == '':
                        clnorigin = '.'
                elif "CLNSRC=" in findInfo[i]:
                    refinfo = findInfo[i]
                elif "CLNDSDB" in findInfo[i] and '.' not in findInfo[i]:
                    refinfo = refinfo +';' + findInfo[i]
                elif "MUT" in findInfo[i]:
                    refinfo = refinfo +';' + findInfo[i]
                
                
            if "CLNSRC=." in refinfo:
                refinfo = refinfo.replace("CLNSRC=.", '.').replace(".;", '')
            
            text_tokens.insert(5, clinsig)
            text_tokens.insert(6, clnbn)
            text_tokens.insert(7, clnorigin)
            text_tokens.insert(8, refinfo)
            
            
            ref = text_tokens[3]
            var = text_tokens[4]
            if len(ref) == 1 and (len(var) == 1 or 'or' in var) or len(ref) <= len(var):
                text_tokens.insert(2, (text_tokens[1]))
            if len(ref) > len(var) and ('or' not in var):
                text_tokens.insert(2, (int(text_tokens[1]) +1))

            text_tokens.remove(text_tokens[3])
            text_tokens = str(text_tokens [0:9])
            #print text_tokens
            text_tokens = text_tokens.replace('[','').replace(']', '').replace("'", "").replace(",", "\t")
            text_tokens = text_tokens.replace(' or ', ',').replace( ".,", '')
            fout.write(text_tokens)
            fout.write('\n')
            
        else: 
            fout.write(line)
            

    print "brca1count=", brca1count, '\n', "brca2count=", brca2count
    fin.close()
    fout.close()
    

VCF_Modification('clinvar_20131203.vcf', 'refined_clinvar.vcf')