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
        out_file.write('chr , start, end, ref, var, clinsig, references, exon/intron, occurrence of mutation')
        out_file.write('\n')
        seen = set() # set for fast O(1) amortized lookup
        
        for row in data:
            impKey = str(row[2:5])
            if impKey in seen: 
                continue # skip duplicate
            seen.add(impKey)
            row = (str(row)).replace(']', '').replace('[','').replace("'",'')
            row = row + ', ' + str(counterDict[impKey])
            out_file.write(row)
            out_file.write('\n')
            


count_and_remove_duplicates('brca1_data.csv', 'brca1_dataB.csv')
count_and_remove_duplicates('brca2_data.csv', 'brca2_dataB.csv')  
