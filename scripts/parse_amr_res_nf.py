#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:16:12 2023

@author: Hakan
NEED TO EDIT TO BE ABLE TO TAKE INPUT FILES FROM COMMAND LINE
AND WORK WITH NEXTFLOW
"""
import os
import pandas as pd
import json
import ast
import argparse

#Define command-line input
parser = argparse.ArgumentParser()
parser.add_argument("-rgi", "--rgi_dir_path", dest = "rgi", help = "Path to directory with RGI results")
parser.add_argument("-abr", "--abricate_file", dest = "abr", help = "Abricate results tab")
parser.add_argument("-chm", "--checkm-file", dest = "chm", help = "Parsed summary of Checkm results")
parser.add_argument("-tax", "--taxonomy-file", dest = "tax", help = "GTDB-Tk summary tsv") 
parser.add_argument("-o", "--output", dest = "output", help = "Name of output file, include suffix")
args = parser.parse_args()

#File and directory paths
index_file = r"/home/jay/data/amr_finding/aro_index.tsv"
rgi_dir_path = args.rgi
abricate_file = args.abr
checkm_file = args.chm
taxonomy_file = args.tax
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
#abricate['#FILE']=abricate['#FILE'].str.replace(r'bins/','')
abricate['#FILE']=abricate['#FILE'].map(os.path.basename)

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

#make a file with bin_gene(bin name + seq number), aro, tool, &identity/Best_Identities (depending on tool)
print('=====Creating aro_info file======')
meta_abr =abricate.loc[:,('SEQUENCE','ARO', '%IDENTITY')]
meta_abr['TOOL'] = 'ABR'
meta_rgi = rgi.loc[:,('SEQUENCE','ARO')]
meta_rgi['%IDENTITY'] = rgi.loc[:,('Best_Identities')]
meta_rgi['TOOL'] = 'RGI'
meta = pd.concat([meta_abr,meta_rgi],ignore_index=True)
meta.to_csv('erken_aro_info.csv', index=False)

print('=====Reading checkm data======')
with open(checkm_file) as f:
    data = f.read()

checkm = ast.literal_eval(data)

print('=====Reading GTDB-Tk data======')
taxonomy = pd.read_csv(taxonomy_file, sep='\t')
taxonomy['user_genome'] = taxonomy['user_genome'].str.replace(r'.fna','')

print('=====Combining Abricate and RGI Strict results======')
#extract file and aro columns, combine to one dataframe
#goal: make a json file
rgi_tofile = rgi[['#FILE','ARO']][rgi['Cut_Off']!='Loose']
abr_tofile=abricate[['#FILE','ARO']]
comb_tofile = pd.concat([abr_tofile,rgi_tofile],ignore_index=True)

print('Creating nested dictionary')
#bin names and aros in dictionary
comb_tofile['#FILE'] = comb_tofile['#FILE'].str.replace(r'.ffn.gz', '')
bins = comb_tofile['#FILE'].unique()
#one dict for all bins
dic = {}
#individual dicts for each bin, save in dic
for b in range(len(bins)):
    print(f'Adding bin: {bins[b]}')
    if taxonomy[taxonomy['user_genome']==bins[b]].empty:
        tax = 'Unclassified'
    else:
        tax = taxonomy['classification'][taxonomy['user_genome']==bins[b]].squeeze()
    t_dic = {'Bin_Id': bins[b], 'ARO': comb_tofile['ARO'][comb_tofile['#FILE']==bins[b]].to_list(),
             'Completeness': checkm[bins[b]]['Completeness'], 'Contamination': checkm[bins[b]]['Contamination'],
             'Lineage': tax}
    dic[bins[b]]=t_dic

print('Saving dict as json file')
with open(args.output, "w") as out:
    json.dump(dic, out, indent=4) 
