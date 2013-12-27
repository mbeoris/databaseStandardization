'''
Created on Dec 26, 2013

@author: mbeoris
'''

#v1.1
#changed the del function to have ref of deleted nt and var of '-'

#input will be a varscan file, if there is only a header 
#skip the file (assuming one line case)
#if -Nt modify ref and end position

import os

#for SNP varscan files, does not include ins/del
def SNP_varscan_modification(inputfile, newfile):
    
    ftest = open(inputfile)
    num_lines = sum(1 for line in ftest)
    ftest.close()
    
    fin = open(inputfile, 'rU')
    #only creates new file if there is more than a header
    if num_lines > 1:
        fout = open(newfile, 'w')
    
    for line in fin:
        if num_lines == 1:
            continue
        elif "Chrom" in line:
            fout.write("chr" + '\t' +  "start" + '\t' + "end"+ '\t' + "ref" + '\t' +"var" \
                       + '\t'+ "Coverage" + '\t'+"AF" + '\t'+ "StrandBias" + '\t' + 'Qual1' \
                       + '\t' + 'Qual2'+ '\t' + 'Pvalue' + '\n')
        else:
            line = line.split('\t')
            strandBias = 1.0*(min (int(line[16]), int(line[17])))/(int(line[16]) + int(line[17]))
            line[0] = str(line[0]).replace('chr', '')
            line[5] = int(line[4])+int(line[5])
            line[4] = str(line[18]).replace('\n', '')
            line.insert(2, (line[1]))
            line.pop(4)
            line.insert(7,strandBias)
            line.pop(8)
            line.pop(8)
            line = str(line[0:11]).replace("'", '').replace( ',', '\t' ).replace('[','').replace(']', '')
            if os.path.exists(newfile):
                fout.write(line)
                fout.write('\n')
    fin.close()
    if os.path.exists(newfile):
        fout.close()
        
              
def indel_varscan_modification(inputfile, newfile):
    
    ftest = open(inputfile)
    num_lines = sum(1 for line in ftest)
    ftest.close()
    
    fin = open(inputfile, 'rU')
    fout = open(newfile, 'w')
    
    for line in fin:
        if num_lines == 1:
            fin.close()
            fout.close()
            continue
        
        elif "Chrom" in line:
            fout.write("chr" + '\t' +  "start" + '\t' + "end"+ '\t' + "ref" + '\t' +"var" \
                       + '\t'+ "Coverage" + '\t'+"AF" + '\t'+ "StrandBias" + '\t' + 'Qual1' \
                       + '\t' + 'Qual2'+ '\t' + 'Pvalue' + '\n')
        
        else:
            line = line.split('\t')
            strandBias = 1.0*(min (int(line[16]), int(line[17])))/(int(line[16]) + int(line[17]))
            line[0] = str(line[0]).replace('chr', '')
            #sum to get depth
            line[5] = int(line[4])+int(line[5])
            line[4] = str(line[18]).replace('\n', '')
            #replicating start site
            line.insert(2, (line[1]))
            line.pop(4)
            line.insert(7,strandBias)
            line.pop(8)
            line.pop(8)
            
            #handling delete cases
            if '-' in line[4]:
                delnt = str(line[4]).replace('-','')
                var = '-'
                ref = delnt
                startpos = int(line[1]) + len(delnt)
                line[1] = line[2] = startpos
                line[3] = ref
                line[4] = var
            
            #handling insertion cases
            elif '+' in line[4]:
                insnt = str(line[4]).replace('+', '')
                var = str(line[3]) + insnt
                line[4] = var
            
            
            line = str(line[0:11]).replace("'", '').replace( ',', '\t' ).replace('[','').replace(']', '')
            fout.write(line)
            fout.write('\n')
                        
            
                
            

    fin.close()
    fout.close()


#SNP_varscan_modification('varscan.SNP.vcf', "refinedvarscan.SNP.vcf")
#indel_varscan_modification('varscan.indel.vcf', 'refinedvarscan.indel.vcf')       




