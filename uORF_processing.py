# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import re
import numpy as np
import pandas as pd

with open("/Users/mingchuxu/Downloads/sequence_only.txt") as f:
    seq_list = f.readlines()
    
for i in range(len(seq_list)):
    
    seq_list[i] = seq_list[i].rstrip()



def search_start(seq):

	a = next((i for i,v in enumerate(seq) if v == 'ATG'), None)

	if a == None:
		return
	else:
		if len(array) == 0:
			array.append(a)
		else:
			array.append(a + array[-1])
		
		search_end(seq[a:])

def search_end(seq):

	b = next((i for i,v in enumerate(seq) if (v == 'TAG') | (v == 'TAA') | (v == 'TGA')), None)

	if b == None:
		array.append(-1)
		return
	else:
		array.append(b + array[-1])
		search_start(seq[b:])

def convert(array, frame):

    a = np.reshape(np.array(array),(-1,2))
    m = a[:,0] * 3 + 1 + frame
    n = a[:,1] * 3 + 3 + frame
    return np.column_stack((m, n))
	


## Working code


mrna_coor_0 = []
mrna_coor_1 = []
mrna_coor_2 = []

for i in range(0, len(seq_list)):
    
    array = []

    tpl = re.findall('...?', seq_list[i][0:])

    search_start(tpl)

    mrna_coor_0.append(convert(array, 0))
    
    
for i in range(0, len(seq_list)):
    
    array = []

    tpl = re.findall('...?', seq_list[i][1:])

    search_start(tpl)

    mrna_coor_1.append(convert(array, 1))
    
    
for i in range(0, len(seq_list)):
    
    array = []

    tpl = re.findall('...?', seq_list[i][2:])

    search_start(tpl)

    mrna_coor_2.append(convert(array, 2))
    

data = pd.read_csv("/Users/mingchuxu/Downloads/without_sequence.txt", sep="\t", header=None)




def coor_to_array(str1, str2):
    
    a1 = str1.split(',')
    b1 = [ int(x) for x in a1 ]
    
    a2 = str2.split(',')
    b2 = [ int(x) for x in a2 ]
     
    return np.column_stack((b1, b2))



exon_coor = []

for i in range(0, data.shape[0]):
    
    exon_coor.append(coor_to_array(data[6][i], data[7][i]))
    

def mrna_to_genome(M, a, sign):
    
    if sign == "+":
        
        k = 1
        
    else:
        
        k = 0
        
        
    if a == 0:
            
        return M[-1][k]
    
    else:
            
        p = M[:,1] - M[:,0] + 1
            
        q = np.array([0, p[0]])
            
        for i in range(1, len(p)):
                
            q = np.append(q, q[i] + p[i])
                
        j = 0
            
        while (q[j] < a):
                
            j = j + 1
            
        diff = a - q[j-1]
            
        return M[j-1][1-k] + (k * 2 - 1) * (diff - 1)
    

def array_to_genome(M, N, sign):
    
    nrow = np.shape(N)[0]
    
    a = np.array([])
    
    a = a.astype(int)
    
    for i in range(nrow):
        
        a = np.append(a, [mrna_to_genome(M, j, sign) for j in N[i]])
        
    return np.reshape(a, (-1,2))




orf_coor_0 = []
orf_coor_1 = []
orf_coor_2 = []

for i in range(len(exon_coor)):
    
    orf_coor_0.append(array_to_genome(exon_coor[i], mrna_coor_0[i], data[5][i]))
    orf_coor_1.append(array_to_genome(exon_coor[i], mrna_coor_1[i], data[5][i]))
    orf_coor_2.append(array_to_genome(exon_coor[i], mrna_coor_2[i], data[5][i]))
    
    
UTR_length = []

for i in range(len(seq_list)):
    
    UTR_length.append(len(seq_list[i]))
    

new_data = pd.DataFrame(
        
        {'UTR_length': UTR_length,
         'exon_coor': exon_coor,
         'mrna_coor_0': mrna_coor_0,
         'mrna_coor_1': mrna_coor_1,
         'mrna_coor_2': mrna_coor_2,
         'orf_coor_0': orf_coor_0,
         'orf_coor_1': orf_coor_1,
         'orf_coor_2': orf_coor_2,
         'seq':seq_list
        })
    
uorf_data = pd.concat([data, new_data], axis=1)

uorf_data.drop(uorf_data.columns[6:8], axis = 1, inplace = True)    

colnames = uorf_data.columns.values

colnames[0:6] = ['ID', 'symbol', 'chr', 'start', 'end', 'strand']    

uorf_data.columns = colnames
    
    
    
        
        

        
        
        
        
        
        
    


    




    


