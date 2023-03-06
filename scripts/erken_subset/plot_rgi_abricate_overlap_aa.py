#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:16:12 2023

@author: Hakan
Messy script to compare RGI and Abricate to CARD
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

#File and directory paths
os.chdir('/home/jay/work_dir/python')
index_file=r"/home/jay/data/erken_results/03_rgi/card_stuff/aro_index.tsv"
rgi_dir_path=r"/home/jay/data/erken_results/03_rgi/res_aa"
abricate_file=r"/home/jay/data/erken_results/02_abricate/card_coding_50results.tab"

"""
a function that takes an accession and gene, then finds which aro it
corresponds to. Also index should be a file with ARO accessions
and ARO names.
returns aro as int
"""
def find_aro(gene, index):
    aro = index.loc[index['ARO Name']==gene, 'ARO Accession'].squeeze()
    mech = index.loc[index['ARO Accession']==aro, 'Resistance Mechanism'].squeeze()
    return int(aro), mech

print('=====Reading and formatting index======')
#format index file
#col_gone = ['CVTERM ID', 'Model Sequence ID','Model ID', 'Model Name', 'Protein Accession',
#            'AMR Gene Family','Drug Class', 'Resistance Mechanism', 'CARD Short Name']
index = pd.read_csv(index_file, sep='\t', header=0)
#index = index.drop(col_gone, axis=1)
index['ARO Name'] = index['ARO Name'].str.replace(r' ', '_')
index['ARO Accession'] = index['ARO Accession'].str.replace(r'ARO:', '')

print('====== Reading and formatting abricate results======')
abricate = pd.read_csv(abricate_file, sep='\t', header=0)
abricate['#FILE']=abricate['#FILE'].str.replace(r'bins/','')

print('Adding ARO terms')
# Add ARO terms to abricate results
#dtype for ARO in rgi is int64
#also adding the resistance mechanism
abricate['ARO']=0
abricate['Mechanism']=''
for ind in abricate.index:
    aro = find_aro(abricate['GENE'][ind], index)
    abricate.loc[ind, 'ARO'] = aro[0]
    abricate.loc[ind, 'Mechanism'] = aro[1]

print('======Reading and formatting RGI results====')
#RGI-results
rgi = pd.DataFrame()
rgi = pd.concat([pd.read_csv(rgi_dir_path+'//'+file, sep='\t', header=0)
          for file in os.listdir(rgi_dir_path) if file.endswith('.txt')],
                ignore_index=True)
#rgi_del_col = ['ORF_ID', 'Model_type','Predicted_DNA','Predicted_Protein',
#               'CARD_Protein_Sequence','Nudged', 'Note']
#rgi = rgi.drop(rgi_del_col, axis=1)
# give rgi a bin name column matching abricates format
rgi['#FILE']=''
rgi['SEQUENCE']=''
for ind in rgi.index:
    bin = rgi['ORF_ID'][ind].split(' ',1)[0]
    bin = bin.rsplit('_',1)[0]
    rgi.loc[ind, '#FILE'] = bin + '.ffn.gz'
    seq = rgi['ORF_ID'][ind].split(' ')
    rgi.loc[ind, 'SEQUENCE']=seq[0]
    print(bin + ': Done')

# Takes two dataframes, plots a venn diagram of specified column, saves as png
def plot_venn(df1, df2, df1_name, df2_name, fig_name, column, title=''):
    if not df1.empty and not df2.empty:
        plt.figure(figsize=(4,4))
        set1 = set(df1[column])
        set2 = set(df2[column])
        venn2([set1, set2], (df1_name, df2_name))
        plt.title(title)
        plt.savefig(fig_name+'.png')
        plt.show()
    else:
        print('One or both of the dataframes are empty.')

# plot venn diagrams
print('=====Plotting first venn=====')
plot_venn(abricate, rgi, 'Abricate', 'RGI', 'venn_abr_rgi_aa','ARO', 'AA ARO overlap')
print('=====Plotting Strict venn======')
plot_venn(abricate, rgi[rgi['Cut_Off']=='Strict'], 'Abricate', 'RGI Strict', 'venn_abr_rgi_strict_aa','ARO', 'Strict AA ARO overlap')
print('======Plotting Perfect venn=====')
plot_venn(abricate, rgi[rgi['Cut_Off']=='Perfect'], 'Abricate', 'RGI Perfect', 'venn_abr_rgi_perfect_aa','ARO', 'Perfect AA ARO overlap')

print('=====Plotting bin venns=====')
plot_venn(abricate, rgi, 'Abricate', 'RGI', 'bin_venn_abr_rgi_aa','#FILE','AA Bin overlap')
print('=====Strict======')
plot_venn(abricate, rgi[rgi['Cut_Off']=='Strict'], 'Abricate', 'RGI Strict', 'bin_venn_abr_rgi_strict_aa','#FILE','Strict AA bin overlap')
print('======Perfect=====')
plot_venn(abricate, rgi[rgi['Cut_Off']=='Perfect'], 'Abricate', 'RGI Perfect', 'bin_venn_abr_rgi_perfect_aa','#FILE', 'Perfect AA bin overlap')

print('=====Plotting Sequence venns=====')
plot_venn(abricate, rgi, 'Abricate', 'RGI', 'seq_venn_abr_rgi_aa','SEQUENCE','AA Sequence overlap')
plot_venn(abricate, rgi[rgi['Cut_Off']=='Strict'], 'Abricate', 'RGI Strict', 'seq_venn_abr_rgi_strict_aa','SEQUENCE','Strict AA Sequence overlap')

print('=====Plotting Best hits venns=====')
plot_venn(abricate[abricate['%IDENTITY']>75], rgi[rgi['Cut_Off']=='Strict'], 'Abricate', 'RGI Strict', 'bin_venn_abr75_rgi_strict_aa','#FILE','AA best bins overlap')
plot_venn(abricate[abricate['%IDENTITY']>75], rgi[rgi['Cut_Off']=='Strict'], 'Abricate', 'RGI Strict', 'aro_venn_abr75_rgi_strict_aa','ARO','AA best ARO overlap')

plot_venn(abricate[abricate['%IDENTITY']>75], rgi[rgi['Cut_Off']=='Strict'], 'Abricate', 'RGI Strict', 'seq_venn_abr75_rgi_strict_aa','SEQUENCE','AA best Sequence overlap')


#Find differences
print('===== Finding differences =====')
aro_abr = abricate['ARO'].to_list()
aro_rgi = rgi['ARO'].to_list()
s = set(aro_rgi)
aro_diff = [x for x in aro_abr if x not in s]

print('AROs that are in abricate res but not in rgi res:')
print(aro_diff)

aro_abr = abricate['ARO'].to_list()
aro_rgi = rgi['ARO'][rgi['Cut_Off']=='Strict'].to_list()
s = set(aro_rgi)
aro_sim = [x for x in aro_abr if x in s]

print('AROs that are in common between abricate and RGI Strict:')
print(aro_sim)

