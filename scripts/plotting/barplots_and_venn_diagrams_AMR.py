# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:14:58 2023

@author: Jay Håkansson
As a script this is very messy, but it should give some insight into
how the figures for my report were made.
"""

import os
import pandas as pd
import numpy as np
from matplotlib_venn import venn2 
from matplotlib import pyplot as plt

# define functions

def format_df(df, dataset_name=''):
    if dataset_name=='':
        raise Exception("You need to set a dataset name")
    df['Sample'] = ''
    df['bin_id'] = ''
    df['ARO_name'] = ''
    df['dataset'] = dataset_name
    for ind in df.index:
        sample = df['Seq_Id'][ind].split('_', 1)[0]
        bin = df['Seq_Id'][ind].rsplit('_', 1)[0]
        df.loc[ind, 'Sample'] = sample
        df.loc[ind, 'bin_id'] = bin
        df.loc[ind, 'ARO_name'] = aro_index.loc[aro_index['ARO Accession']
                                                 == str(df.loc[ind, 'ARO']),
                                                 'CARD Short Name'].item()
    return df

def add_drug_and_resistance(df):
    df['Drug_class']=''
    df['Resistance_Mechanism']=''
    for ind in df.index:
        aro = str(df['ARO'][ind])
        df.loc[ind, 'Drug_class'] = aro_index.loc[aro_index['ARO Accession']
                                                  ==aro, 'Drug Class'].item()
        df.loc[ind, 'Resistance_Mechanism'] = aro_index.loc[aro_index['ARO Accession']
                                                  ==aro, 'Resistance Mechanism'].item()
    return df

def find_named_rank(tax_ranks, lin, t):
    if t==0:
        return lin[0]
    elif len(lin[t]) <= 3:
        return find_named_rank(tax_ranks, lin, t-1)
    else:
        return lin[t]


def format_tax_ranks(df):
    tax_ranks = ['Domain','Phylum','Class','Order','Family','Genus','Species']
    #create columns
    for t in tax_ranks:
        df[t]=''
    #fill columns
    for ind in df.index:
        lin = df['Lineage'][ind].split(';')
        for l in range(0,len(lin)):
            df.loc[ind, tax_ranks[l]] = lin[l]
        if len(lin) < len(tax_ranks):
            for i in range(len(lin), len(tax_ranks)):
                #df.loc[ind, tax_ranks[i]] = 'Unclassified'
                df.loc[ind, tax_ranks[i]] = lin[0]
        elif len(lin[6]) <= 3:
            for t in range(1, len(tax_ranks)):
                if len(lin[t]) <= 3:
                    named_rank = find_named_rank(tax_ranks, lin, t)
                    df.loc[ind, tax_ranks[t]] = named_rank + ';' + lin[t] + 'Unclassified'
    return df

def add_bin_quality(df):
    df['bin_quality']=''
    for ind in df.index:
        if (df['Completeness'][ind] >=90 and df['Contamination'][ind] <=5):
            df.loc[ind, 'bin_quality']='Good'
        elif (df['Completeness'][ind]>=50 and df['Contamination'][ind]<=10):
            df.loc[ind, 'bin_quality']='Medium'
        elif (df['Completeness'][ind]<50 and df['Completeness'][ind]>0):
            df.loc[ind, 'bin_quality']='Bad'
        else:
            df.loc[ind, 'bin_quality']='Rest'
    return df

                
#
def plot_aro_counts_x(df, topnr = 10, title = '', log10 = False, use_names=False):
    if use_names == True:
        col = 'ARO_name'
    else:
        col = 'ARO'
    if log10 == True:
        aro_counts = np.log10(df[col].value_counts())
    else:
        aro_counts = df[col].explode().value_counts()
        
    top_x_aro_counts = aro_counts.copy(deep=True)
    top_x_aro_counts = top_x_aro_counts.drop(top_x_aro_counts.index[topnr:])
    top_x_aro_counts.plot(kind='barh', title = title)

#not done also not working haha    
def plot_aros_by_species_count(df, topnr = 10, title = ''):
    group = df.groupby(['ARO_name', 'Lineage']).size().groupby(level=0).size()
    group.sort_values(ascending=False)[:topnr].plot(kind='barh', title = title)
    return group.sort_values(ascending=False)[:topnr]
    
def plot_lineage_counts(df, topnr = 10, lineage = 'Lineage', title = '', log10 = False, col = 'tab:blue'): 
    if log10 == True:
        lineage_count = np.log10(df[lineage].value_counts())
    else:
        lineage_count = df[lineage].value_counts()
    top_x_lineages = lineage_count.drop(lineage_count.index[topnr:])
    top_x_lineages.plot(kind='barh', title = title, color = col, figsize=(15,15), fontsize=32)
    return top_x_lineages
    
def bin_count_stats(df):
    total_bins = len(df['bin_id'].unique())
    uniq_bins= len(df.drop_duplicates('bin_id', keep=False)) #this drops *all* rows with duplicate bins
    dup_bins=df[df.duplicated('bin_id', keep=False)] #all duplicates = True, which selects all dups
    dup_bins_uniq = len(dup_bins.drop_duplicates('bin_id')) #keeps only one of each bin
    return total_bins, uniq_bins, dup_bins_uniq
 
def pie_bin_quality(df, title=''):
    col_pal=["palegreen", "g", "c", "dodgerblue"]
    #col_pal=["palegreen","dodgerblue","c","g"]
    uniq_bins = df.drop_duplicates('bin_id', keep='first')
    uniq_bins['bin_quality'].value_counts().plot.pie(autopct='%.1f %%', ylabel='',
                                                     title=title,
                                                     colors=col_pal, figsize=(6,6), fontsize=12)

    

#%%
#read in dfs
print('Reading in dfs')
stratfreshDB_res_path=r"D:\Plugg\Master thesis\AMR results with archaea\Stratfresh"
sites_res_path=r"D:\Plugg\Master thesis\AMR results with archaea\SITES"
motus_path = r"D:\Plugg\Master thesis\downloads\buck_master_table.csv"
aro_index_file = r"D:\Plugg\Master thesis\downloads\amr_finding\aro_index.tsv"

stratfreshDB = pd.DataFrame()
stratfreshDB = pd.concat([pd.read_json(stratfreshDB_res_path+'\\'+file, orient = 'index')
                for file in os.listdir(stratfreshDB_res_path) if file.endswith('.json')],
                ignore_index=True)
stratfreshDB = stratfreshDB.explode('ARO', ignore_index=True)
print('Stratfresh done')

sites = pd.DataFrame()
sites = pd.concat([pd.read_json(sites_res_path+'\\'+file, orient = 'index')
                for file in os.listdir(sites_res_path) if file.endswith('.json')],
                ignore_index=True)
sites = sites.explode('ARO', ignore_index=True)
print('SITES done')

#read in ARO info file and format to match result naming
aro_index = pd.read_csv(aro_index_file, sep='\t')
aro_index['ARO Name'] = aro_index['ARO Name'].str.replace(r' ', '_')
aro_index['ARO Accession'] = aro_index['ARO Accession'].str.replace(r'ARO:', '')

sites_meta=pd.read_csv(r"D:\Plugg\Master thesis\SITES\sites_meta.csv", sep=',')
#they have dumb names. Need to fix to be compatible. Easiest manually tbh.
#sites_meta['sample_name']=sites_meta['label']
#sites_meta.to_csv(r"D:\Plugg\Master thesis\SITES\sites_meta.csv", sep=',', index=False)


#%%
print('Format dfs. This takes a while')
stratfreshDB = format_df(stratfreshDB, 'StratfreshDB')
sites = format_df(sites, 'SITES')
print('Initial formatting done, starting taxonomy formatting')
stratfreshDB = format_tax_ranks(stratfreshDB)
sites = format_tax_ranks(sites)
print('Add bin quality')
stratfreshDB = add_bin_quality(stratfreshDB)
sites = add_bin_quality(sites)


#all that arent single samples are lake coassemblies for this dataset
print('Add sample type')
sites['sample_type']='lake_coasm'
for ind in sites.index:
    #if starts with Sample = single assembly
    if sites['Sample'][ind].startswith('Sample-'):
        sites.loc[ind, 'sample_type']= sites_meta[sites_meta['sample_name']==
                                                  sites['Sample'][ind]]['type'].item()
        
sites= add_drug_and_resistance(sites)
stratfreshDB = add_drug_and_resistance(stratfreshDB)

#%% venn diagram of drug classes
sf_drugs = set(stratfreshDB['Drug_class'].unique())
sites_drugs= set(sites['Drug_class'].unique())
sf_drugs-sites_drugs
sites_drugs-sf_drugs
overlap = sf_drugs&sites_drugs
len(overlap)

A = sf_drugs
B = sites_drugs
v = venn2([A, B], ('StratfreshDB families', 'SITES families'),
          set_colors=('dodgerblue','green'))
v.get_label_by_id('10').set_text('\n'.join(A-B))
v.get_label_by_id('11').set_text('\n'.join(A&B))
v.get_label_by_id('01').set_text('\n'.join(B-A))
plt.show()



#%% venn diagram of amr gene overlap
# sites_amr=plot_aros_by_species_count(sites, topnr = 20, title = 'sites aro by species')
# stratfresh_amr=plot_aros_by_species_count(stratfreshDB, topnr = 20, title = 'stratfresh aro by species')
# sites_set=set(pd.Series.tolist(sites_amr.index))
# stratfresh_set=set(pd.Series.tolist(stratfresh_amr.index))

# A = stratfresh_set
# B = sites_set

# overlap = '\n'.join(A&B)
# overlap= overlap.split('\n')
# overlap.sort()

# v = venn2([A, B], ('StratfreshDB AMR genes', 'SITES AMR genes'),
#           set_colors=('dodgerblue','green'))
# plt.annotate('\n'.join(overlap), xy=v.get_label_by_id('11').get_position() +
#               np.array([0, 0.2]), xytext=(65,-145), ha='center',
#               textcoords='offset points',
#               bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.1))
# plt.show()





#%% venn diagrams of family overlap
# print('Only count each bin once for unique dfs')
# uniq_sites = sites.drop_duplicates('bin_id', keep='first').copy()
# uniq_stratfreshDB = stratfreshDB.drop_duplicates('bin_id', keep='first').copy()

# print('Create the family+class thing')
# uniq_sites['Class_Family']=  uniq_sites['Class'] + ';' + uniq_sites['Family']
# top_sites_familys=plot_lineage_counts(uniq_sites, topnr=20, lineage='Class_Family', col = 'darkseagreen')
# uniq_stratfreshDB['Class_Family']=  uniq_stratfreshDB['Class'] + ';' + uniq_stratfreshDB['Family']
# top_stratfresh_familys=plot_lineage_counts(uniq_stratfreshDB, topnr=20, lineage='Class_Family', col = 'powderblue')

# sites_set=set(pd.Series.tolist(top_sites_familys.index))
# stratfresh_set=set(pd.Series.tolist(top_stratfresh_familys.index))


# venn2([stratfresh_set, sites_set], ('StratfreshDB families', 'SITES families'),
#          set_colors=('dodgerblue','green'))
# plt.title('Overlap between 20 most common families')

# plt.show()

# v= venn2([stratfresh_set, sites_set], ('StratfreshDB families', 'SITES families'),
#          set_colors=('dodgerblue','green'))

# A = stratfresh_set
# B = sites_set
# v = venn2([A, B], ('StratfreshDB families', 'SITES families'),
#          set_colors=('dodgerblue','green'))
# v.get_label_by_id('10').set_text('\n'.join(A-B))
# v.get_label_by_id('11').set_text('\n'.join(A&B))
# v.get_label_by_id('01').set_text('\n'.join(B-A))
# plt.show()

# A = stratfresh_set
# B = sites_set
# v = venn2([A, B], ('StratfreshDB families', 'SITES families'),
#          set_colors=('dodgerblue','green'))
# plt.annotate(',\n'.join(A&B), xy=v.get_label_by_id('11').get_position() +
#              np.array([0, 0.2]), xytext=(-20,40), ha='center',
#              textcoords='offset points',
#              bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
#              arrowprops=dict(arrowstyle='->',              
#                              connectionstyle='arc',color='gray'))
# plt.show()

# overlap = '\n'.join(A&B)
# overlap= overlap.split('\n')
# overlap.sort()

# v = venn2([A, B], ('StratfreshDB families', 'SITES families'),
#          set_colors=('dodgerblue','green'))
# plt.annotate('\n'.join(overlap), xy=v.get_label_by_id('11').get_position() +
#              np.array([0, 0.2]), xytext=(-0,-30), ha='center',
#              textcoords='offset points',
#              bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.1))
# plt.show()



