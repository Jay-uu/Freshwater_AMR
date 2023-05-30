import os
import pandas as pd
import json
import ast

stratfresh_master_path="/home/moritz/proj_folder/uppstore2018116/webexport/stratfreshdb/metadata/master_table.csv"
sites_checkm="/crex/proj/fume/private/jay/process_sites/08_amr_finding/01_checkm/"

print("===Reading stratfresh-file===")
stratfresh=pd.read_csv(stratfresh_master_path, sep=',', header=0)
sf = stratfresh[['bin_id', 'completeness', 'contamination']]
sf.columns= ['Bin Id', 'Completeness', 'Contamination']

print("===Reading sites checkm jsons")
sites=pd.DataFrame()
for file in os.listdir(sites_checkm):
    if file.endswith('.json'):
        with open(sites_checkm + '/' +file) as f:
            data = f.read()
        checkm = ast.literal_eval(data)
        df = pd.DataFrame.from_dict(checkm, orient='index')
        sites = pd.concat([sites, df])


sites = sites[['Bin Id', 'Completeness', 'Contamination']]
both = pd.concat([sites, sf])

both.to_csv("sites_stratfresh_checkm_short.tsv", sep='\t', index = False)
