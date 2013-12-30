'''
Created on Dec 19, 2013

@author: mbeoris
'''

import csv
from collections import Counter

def count_and_remove_duplicates(infile, outfile):
    out = open(infile)
    data = csv.reader(out)
    data.next()
    data=[row for row in data]
    out.close()
    formations=[]
    
    for row in data:
        if row[2]== '':
            continue
        else:
            formations.append(str(row[2:5]))
    
    counterDict = Counter(formations)
    
    
    with open(outfile,'w') as out_file:
        out_file.write('chr \t start \t end\t ref\t var\t clinsig\t references\t exon/intron \t occurrence of mutation')
        out_file.write('\n')
        seen = set() # set for fast O(1) amortized lookup
        
        for row in data:
            impKey = str(row[2:5])
            if impKey in seen: 
                continue # skip duplicate
            seen.add(impKey)
            row = (str(row)).replace(']', '').replace('[','').replace("'",'').replace(',', '\t')
            row = row + '\t ' + str(counterDict[impKey])
            out_file.write(row)
            out_file.write('\n')
            


count_and_remove_duplicates('brca1_data.csv', 'brca1_dataB.vcf')
count_and_remove_duplicates('brca2_data.csv', 'brca2_dataB.vcf')  
