'''
Created on Dec 12, 2013

@author: mbeoris
'''
#goal of this module is to find entries from clinvar_YYYYMMDD.vcf file
#input values will be clinvar filename, the newfile will be a condensed version of
#the original clinvar file


def ClinvarVCF_Modification(inputfile,newfile):
    clinsig =''
    clnbn =''
    clnorigin =''
    refinfo =''
    fin = open(inputfile, 'rU')
    fout = open(newfile, 'w')

    for line in fin:
        if '#CHROM' in line:
            fout.write("chr" + '\t' +  "start" + '\t' + "end"+ '\t' + "ref" + '\t' +"var"+ '\t'+ "CLINSIG" + '\t'+"CLNBN" + '\t'+ "CLNORIGIN"+ '\t' +"INFO" + '\t' + "Multiple Alleles")
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
                text_tokens.insert(2, (int(text_tokens[1]) +len(ref)-1))
        
            
            text_tokens.remove(text_tokens[3])
            

            if 'or ' in var:
                line1 = text_tokens[0:9]
                line2 = text_tokens[0:9]
                for i in range(0, len(line1)):
                    if ' or '  in text_tokens[i]:
                        tmp = str(text_tokens[i])
                        countOr= tmp.count('or')
                        orIndex = i
                        if countOr<2:
                            for i in range(0, len(tmp)):
                                if tmp[i:i+2] == 'or':
                                    orPos = i
                            var1 = tmp[0:orPos]                    
                            var2 = tmp[orPos+2:]
                            line1[orIndex] = var1.replace(" ",'')
                            line2[orIndex] = var2.replace(" ",'')
                        if countOr>1:
                            multOrs = tmp.split(';')
                            multOr1 = ""
                            multOr2 = ""
                            for i in range(0, len(multOrs)):
                                if ' or ' in multOrs[i]:
                                    subUnit = multOrs[i]
                                    for i in range(0, len(subUnit)):
                                        if subUnit[i] == '=':
                                            identifier = subUnit[0:i+1]
                                        if subUnit [i:i+2] == 'or':
                                            orPos = i
                                    var1 = subUnit[0:orPos]
                                    var2 = identifier + subUnit[orPos+2:]
                                multOr1 += ";" + var1 
                                multOr2 += ";" + var2
                            multOr1 = multOr1[1:].replace(" ",'')
                            multOr2 = multOr2[1:].replace(" ",'')
                            line1[orIndex] = multOr1
                            line2[orIndex] = multOr2
                
                  
                line2.append("1")
                line2 = str(line2).replace('[','').replace(']', '').replace("'", "").replace(",", "\t")                                            
                fout.write(line2)                   
                fout.write('\n')                                                                 


            else:
                text_tokens[9] = "0"
                text_tokens = str(text_tokens [0:10])
                text_tokens = text_tokens.replace('[','').replace(']', '').replace("'", "").replace(",", "\t")
                text_tokens = text_tokens.replace(' or ', ',').replace( ".,", '')
                fout.write(text_tokens)
                fout.write('\n')
            
        else: 
            fout.write(line)
            

    fin.close()
    fout.close()
    

ClinvarVCF_Modification('clinvar_20131203raw.vcf', 'refined_clinvar.vcf')
