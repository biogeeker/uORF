#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 13:23:34 2018

@author: mingchuxu
"""

import re
import numpy as np
import pandas as pd
import mygene

with open("/Users/mingchuxu/Downloads/sequence_only.txt") as f:
    seq_list = f.readlines()
    
for i in range(len(seq_list)):
    
    seq_list[i] = seq_list[i].rstrip()
    
    
def find_match(lst, items):
    
    result = [i for i in range(len(lst)) if lst[i] in items]
    
    return result
    

def generate_start_position(x, frame):
    
    tpl = re.findall('...?', x[frame:])
    
    a = find_match(tpl, ["ATG"])
    
    return a

def generate_stop_position(x, frame):
    
    tpl = re.findall('...?', x[frame:])
    
    b = find_match(tpl, ["TAG", "TAA", "TGA"])
    
    return b

def find_next_stop(start, lst):
    
    if len(lst) == 0:
        
        lst.append(-2)
    
    if start > lst[-1]:
    
        return [start, -2]
    
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
    
    start_index = generate_start_position(seq, frame)
    
    stop_index = generate_stop_position(seq, frame)
    
    aa_position = generate_orf_position(start_index, stop_index)
    
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
        
        
    if a <= 0:
            
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
    
    

mrna_id_list = data[0].tolist() 

mg = mygene.MyGeneInfo()   

id_data = mg.querymany(mrna_id_list, scopes="refseq", fields=["uniprot"], 
                 species="human", as_dataframe=True)

    
    
new_data = pd.DataFrame(
        
        {
         'gene_id': id_data['_id'].tolist(),       
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




    

bed_temp = []

for i in range(len(seq_list)):
    
    if len(uorf_data['orf_coor_0'][i]) > 0:
        
        for j in range(len(uorf_data['orf_coor_0'][i])):
            
            dict1 = {}
            
            x = 1 if uorf_data['strand'][i] == "+" else 0
            
            dict1.update({0: uorf_data['chr'][i], 
                          1: uorf_data['orf_coor_0'][i][j][1-x],
                          2: uorf_data['orf_coor_0'][i][j][x],
                          3: uorf_data['strand'][i],
                          4: uorf_data['symbol'][i],
                          5: uorf_data['gene_id'][i],
                          6: uorf_data['ID'][i]})
    
            bed_temp.append(dict1)
            
    if len(uorf_data['orf_coor_1'][i]) > 0:
        
        for j in range(len(uorf_data['orf_coor_1'][i])):
            
            dict1 = {}
            
            x = 1 if uorf_data['strand'][i] == "+" else 0
            
            dict1.update({0: uorf_data['chr'][i], 
                          1: uorf_data['orf_coor_1'][i][j][1-x],
                          2: uorf_data['orf_coor_1'][i][j][x],
                          3: uorf_data['strand'][i],
                          4: uorf_data['symbol'][i],
                          5: uorf_data['gene_id'][i],
                          6: uorf_data['ID'][i]})
            
            bed_temp.append(dict1)
            
    if len(uorf_data['orf_coor_2'][i]) > 0:
        
        for j in range(len(uorf_data['orf_coor_2'][i])):
            
            dict1 = {}
            
            x = 1 if uorf_data['strand'][i] == "+" else 0
            
            dict1.update({0: uorf_data['chr'][i], 
                          1: uorf_data['orf_coor_2'][i][j][1-x],
                          2: uorf_data['orf_coor_2'][i][j][x],
                          3: uorf_data['strand'][i],
                          4: uorf_data['symbol'][i],
                          5: uorf_data['gene_id'][i],
                          6: uorf_data['ID'][i]})
            
            bed_temp.append(dict1)
            

            
exon_temp = []

for i in range(len(seq_list)):
        
        for j in range(len(uorf_data['exon_coor'][i])):
            
            dict1 = {}
                
            dict1.update({0: uorf_data['chr'][i], 
                          1: uorf_data['exon_coor'][i][j][0],
                          2: uorf_data['exon_coor'][i][j][1],
                          3: uorf_data['strand'][i],
                          4: uorf_data['symbol'][i],
                          5: uorf_data['gene_id'][i],
                          6: uorf_data['ID'][i]})
                
            exon_temp.append(dict1)
            
uorf_bound = pd.DataFrame(bed_temp)

UTR_bed = pd.DataFrame(exon_temp)


uorf_bound = uorf_bound[uorf_bound[1] < uorf_bound[2]]

UTR_bed = UTR_bed[UTR_bed[1] < UTR_bed[2]]

uorf_bound.to_csv("/Users/mingchuxu/Downloads/uorf_bound.bed", 
                  sep = '\t', header = False, index = False)

UTR_bed.to_csv("/Users/mingchuxu/Downloads/UTR.bed", 
                  sep = '\t', header = False, index = False)
    

            
            
            
        


    
    
