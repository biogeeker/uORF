#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 13:23:34 2018

@author: mingchuxu
"""

import re
import numpy as np
import pandas as pd

with open("/Users/mingchuxu/Downloads/sequence_only.txt") as f:
    seq_list = f.readlines()
    
for i in range(len(seq_list)):
    
    seq_list[i] = seq_list[i].rstrip()
    
    
def find_match(lst, items):
    
    result = [i for i in range(len(lst)) if lst[i] in items]
    
    return result
    

def generate_st_position(x, frame):
    
    tpl = re.findall('...?', x[frame:])
    
    a = find_match(tpl, ["ATG"])
    
    b = find_match(tpl, ["TAG", "TAA", "TGA"])
    
    return [a, b]

def find_next_stop(start, lst):
    
    if len(lst) == 0:
        
        lst.append(-1)
    
    if start > lst[-1]:
    
        return [start, -1]
    
    else:
        
        i = 0
        
        while start > lst[i]:
            
            i = i + 1
            
        return [start, lst[i]]


def generate_orf_position(start, stop):
    
    arr = []
    
    if len(start) == 0:
        
        return arr
    
    else:
        
        for i in range(len(start)):
            
            arr.append(find_next_stop(start[i], stop))
            
        return arr
    
def aa_to_cdna(array, frame):
    
    if len(array) == 0:
        
        return []

    else:
    
        a = np.array(array)
        m = a[:,0] * 3 + 1 + frame
        n = a[:,1] * 3 + 3 + frame
        
        return np.column_stack((m, n))

def predict_orf(seq, frame):
    
    st_index = generate_st_position(seq, frame)
    
    aa_position = generate_orf_position(st_index[0], st_index[1])
    
    return aa_to_cdna(aa_position, frame)


mrna_coor_0 = []
mrna_coor_1 = []
mrna_coor_2 = []


for i in range(len(seq_list)):
    
    mrna_coor_0.append(predict_orf(seq_list[i], 0))
    mrna_coor_1.append(predict_orf(seq_list[i], 1))
    mrna_coor_2.append(predict_orf(seq_list[i], 2))
    




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
    
    
new_data = pd.DataFrame(
        
        {'exon_coor': exon_coor,
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


def generate_interval(chrom, coor, strand, symbol, ID):
    
    if len(coor) == 0:
        
        return
    
    else:
        
        DF = pd.DataFrame()
    
        for i in range(len(coor)):
            
            rec = pd.Dataframe([chrom, coor[i][0], coor[i][1], strand, symbol, ID])
            
            DF.append(rec)
        
        return DF
    

bed_file = pd.DataFrame()

for i in range(len(seq_list)):
    
    if len(uorf['orf_coor_0'][i]) > 0:
    
        bed_file.append(generate_interval(uorf['chr'][i], uorf['orf_coor_0'][i]))
    

            
            
            
        


    
    
