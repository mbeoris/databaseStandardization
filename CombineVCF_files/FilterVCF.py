'''
Created on Dec 27, 2013

@author: mbeoris
'''

import os

#filter out bad entries, assuming: low depth, low allele frequency,
#and high Strand Bias are undesirable

def Filter_VCF(inputfile, depthCutOff=0, afCutOff=0, sbCutOff=1):
    wd = os.getcwd()
    outfile = wd + "\\output\\INDELs\\filtered_" + inputfile
    inputfile = wd + '\\output\\INDELs\\' + inputfile
    print inputfile
    print outfile
    if os.path.exists(inputfile):
        fin = open(inputfile,'rU')
        fnew = open(outfile, 'w')
        cov = 5
        sb = 7
        af = 6
        for line in fin:
            splitLine = line.split('\t')
            print (float(str(splitLine[af]).replace("%", '')))
            if (float(splitLine[cov]) < depthCutOff) or (float(splitLine[sb]) > sbCutOff) or ((float(str(splitLine[af]).replace("%", ''))) < afCutOff):
                continue
            else:
                fnew.write(line)
    fin.close()
    fnew.close()
    

Filter_VCF('combined_INDELs.vcf', 100, sbCutOff=.01)
                