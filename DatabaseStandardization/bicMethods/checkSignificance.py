'''
Created on Dec 17, 2013

@author: mbeoris
'''

#need to determine if any intron mutations have clinical significance

def checkIntronSignificance(inputfile, clinSigCol):
    countNumSig =0
    countIntronMut =0
    sigIntrons = []
    fin = open(inputfile, 'rU')
    
    for line in fin:
        text_tokens = line.split('\t')
        if "I" in text_tokens[1]:
            countIntronMut += 1
            print text_tokens[clinSigCol]
            if str(text_tokens[clinSigCol]).upper() == "YES" or str(text_tokens[clinSigCol]).upper()== "UNKNOWN":
                print "FOUND SIGNIFICANCE"
                sigIntrons.append(text_tokens[2])
                countNumSig +=1
    
    print "Intron Mutations: " , countIntronMut
    print countNumSig
    
    print sigIntrons
    fin.close()
    
    

checkIntronSignificance("brca1_data.txt", 13)
    
    
    