'''
Created on Dec 2, 2013

@author: mbeoris

original exctracted information from:
http://brca.iarc.fr/LOVD/variants.php?select_db=BRCA1&action=view_unique&order=Variant%2FDNA%2CASC&hide_col=&show_col=&limit=250

script translates LOVD_BRCA1 from raw text to .xls
prev format: Replaces current DNA in (c.1234T>C) format with just a number,
*after inspection of reference sequence (gb U14680), realized that
correct position came from BIC change column, not DNA change column
then creates two new columns labeled "Old nt" and "New nt"

'''
print "modifyRawText v.1"
#inputFile raw text from site, outputExcel desired name for excel file

def modifyRawText(inputFile, outputExcel, posCol):
    print posCol
    
    rawText = open(inputFile, 'r')
    newExcel = open(outputExcel, 'w')
    DNAchangeIndex = 0;
       
    for line in rawText:
        line = line.replace(',', '.')
        text_tokens = line.split('\t')
        for i in range(0,len(text_tokens)):
            if posCol in text_tokens[i]:
                print i
                DNAchangeIndex = i;
                text_tokens.insert(i+1, "ref")
                text_tokens.insert(i+2, "var")

        if DNAchangeIndex != 0 and text_tokens[DNAchangeIndex] != posCol:
            oldForm = text_tokens[DNAchangeIndex]
            oldNt = oldForm[-3]
            text_tokens.insert(DNAchangeIndex + 1, oldNt)
            newNt = oldForm[-1]
            text_tokens.insert(DNAchangeIndex + 2, newNt)
            if "c." in text_tokens[DNAchangeIndex]:
                text_tokens[DNAchangeIndex] = oldForm[2:-3]
            else:
                text_tokens[DNAchangeIndex] = oldForm[0:-3]
        
        text_tokens = str(text_tokens)
        newLine = text_tokens.replace("'", "").replace('[','').replace(']', '')
        newLine = newLine + ' \n'
        newExcel.write(newLine)
  

    rawText.close()
    newExcel.close()


modifyRawText("LOVD_BRCA1_12.2.13.txt", "LOVD_BRCA1_12.2.13.csv", "BIC DNA change")

modifyRawText("LOVD_BRCA2_12.10.13.txt", "LOVD_BRCA2_12.10.13.csv", "BIC DNA change")

