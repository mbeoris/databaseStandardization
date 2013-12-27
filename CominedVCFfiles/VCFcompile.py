'''
Created on Dec 26, 2013

@author: mbeoris
'''
#v1

#need to access files using IDs in target.txt,
#open target files, and write all information
#into one large file



import os
from varscanModification import SNP_varscan_modification
from varscanModification import indel_varscan_modification
#from FilerVCF import Filter_VCF

#the parentFolder input for all methods will be the 
#folder where the desired target.txt file exists

def get_targets(parentFolder):
    os.chdir(parentFolder)
    targ = open('targets.txt', 'rU')
    sampleNames = set()
    for line in targ:
        line = line.split('\t')
        sampleNames.add(line[1])
    targ.close() 
    return sampleNames



#SNP functions
def modify_all_varSNPs(parentFolder):
    os.chdir(parentFolder)
    wd = os.getcwd()
    sampleNames = get_targets(wd)
    
    for i in sampleNames:
        targPath = wd + '\output\\' + i +'\\varscan.SNP.vcf'
        if os.path.exists(targPath):
            newPath = wd + '\output\\' + i +'\\updated_varscan.SNP.vcf'
            SNP_varscan_modification(targPath, newPath)
            
    
def VCF_SNPcombine(parentFolder):
    #change current directory to the parentFolder path
    os.chdir(parentFolder)
    
    #create new directory if it doesn't already exist
    if not os.path.exists('output/SNPs'): 
        os.makedirs('output/SNPs')
    
    #save working directory, which now should be the same as parentFolder    
    wd = os.getcwd()

    newfile = wd +'\output\SNPs\combined_SNPs.vcf'
    fnew = open(newfile, 'w')
    sampleNames = get_targets(wd)
    for i in sampleNames:
        targPath = wd + '\output\\' + i +'\\updated_varscan.SNP.vcf'
        if os.path.exists(targPath):
            fopen = open(targPath, 'rU')
            for line in fopen:
                if "Pvalue" in line:
                    continue
                else:
                    line = str(line).replace('\n','\t') + i + '\n'
                    fnew.write(line)
            fopen.close() 
    fnew.close()
    
    

#INDEL FUNCTIONS
def modify_all_varINDELs(parentFolder):
    os.chdir(parentFolder)
    wd = os.getcwd()
    sampleNames = get_targets(wd)
    
    for i in sampleNames:
        targPath = wd + '\output\\' + i +'\\varscan.indel.vcf'
        if os.path.exists(targPath):
            newPath = wd + '\output\\' + i +'\\updated_varscan.indel.vcf'
            indel_varscan_modification(targPath, newPath)


def VCF_INDELcombine(parentFolder):
    os.chdir(parentFolder)
    
    if not os.path.exists('output/INDELs'): 
        os.makedirs('output/INDELs')
        
    wd = os.getcwd()

    newfile = wd +'\output\INDELs\combined_INDELs.vcf'
    fnew = open(newfile, 'w')
    sampleNames = get_targets(wd)
    for i in sampleNames:
        targPath = wd + '\output\\' + i +'\\updated_varscan.indel.vcf'
        if os.path.exists(targPath):
            fopen = open(targPath, 'rU')
            for line in fopen:
                if "Pvalue" in line:
                    continue
                else:
                    line = str(line).replace('\n','\t') + i + '\n'
                    fnew.write(line)
            fopen.close()   
    fnew.close()  


wd = os.getcwd()

#modify_all_varSNPs(wd)
#VCF_SNPcombine(wd)
#modify_all_varINDELs(wd)
#VCF_INDELcombine(wd)
