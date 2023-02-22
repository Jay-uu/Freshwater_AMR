#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:16:12 2023

@author: Hakan
Bad messy script to compare RGI and Abricate to CARD
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

#File and directory paths
os.chdir('/home/jay/work_dir/python')
index_file=r"/home/jay/data/erken_results/03_rgi/card_stuff/aro_index.tsv"
rgi_dir_path=r"/home/jay/data/erken_results/03_rgi/res"
abricate_file=r"/home/jay/data/erken_results/02_abricate/card_50results.tab"


"""
a function that takes an accession and gene, then finds which aro it
corresponds to. Also index should be a file with ARO accessions
and ARO names.
returns aro as int
"""
def find_aro(gene, index):
    aro = index.loc[index['ARO Name']==gene, 'ARO Accession'].squeeze()
    return int(aro)

#format index file
col_gone = ['CVTERM ID', 'Model Sequence ID','Model ID', 'Model Name', 'Protein Accession',
            'AMR Gene Family','Drug Class', 'Resistance Mechanism', 'CARD Short Name']
index = pd.read_csv(index_file, sep='\t', header=0)
index = index.drop(col_gone, axis=1)
index['ARO Name'] = index['ARO Name'].str.replace(r' ', '_')
index['ARO Accession'] = index['ARO Accession'].str.replace(r'ARO:', '')

abricate = pd.read_csv(abricate_file, sep='\t', header=0)
abricate['#FILE']=abricate['#FILE'].str.replace(r'bins/','')

# Add ARO terms to abricate results
#dtype for ARO in rgi is int64
abricate['ARO']=0
for ind in abricate.index:
    aro = find_aro(abricate['GENE'][ind], index)
    abricate.loc[ind, 'ARO'] = aro

#RGI-results
rgi = pd.DataFrame()
rgi = pd.concat([pd.read_csv(rgi_dir_path+'//'+file, sep='\t', header=0)
          for file in os.listdir(rgi_dir_path) if file.endswith('.txt')],
                ignore_index=True)
rgi_del_col = ['ORF_ID', 'Model_type','Predicted_DNA','Predicted_Protein',
               'CARD_Protein_Sequence','Nudged', 'Note']
rgi = rgi.drop(rgi_del_col, axis=1)

# Plot venn diagram and save

plt.figure(figsize=(4,4))
set1 = set(abricate['ARO'])
set2 = set(rgi['ARO'])

venn2([set1, set2], ('Abricate', 'RGI'))
plt.savefig('venn_abr_rgi.png')
plt.show()
