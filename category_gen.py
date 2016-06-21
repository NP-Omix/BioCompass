from sys import argv
from Bio import SeqIO
import pandas as pd
import re

script, strain_name = argv

table1_df = []

table1_df = pd.read_csv('%s_table1.csv' % strain_name, sep='\t')

table1_df['product'].fillna('None', inplace=True)

subcluster_df = pd.read_csv('subcluster_dictionary.csv', sep='\t')

subcluster_dict = subcluster_df.set_index('product')['category'].to_dict()

col7 = []
col8 = []
missing_from_dict = []

table1_df['product'] = table1_df['product'].astype(str)

for line in table1_df['product']:
    for key in subcluster_dict:
        m = re.search(key, line, re.I)
        if line != None and m != None:
            col7.append(subcluster_dict[key])
            col8.append(line)

frames = {'category':col7,'product':col8}

new_cols_df = pd.DataFrame(frames, index=None)

table1_df = pd.merge(table1_df, new_cols_df, on='product', how='outer')

table1_df['category'].fillna('hypothetical', inplace=True)

table1_df = table1_df.sort_values('locus_tag',axis=0, ascending=True, inplace=False)

table1_df.drop_duplicates(subset='locus_tag', inplace=True)

table1_df.reset_index(drop=True, inplace=True)

table1_handle = open('%s_table1.csv' % strain_name, "w")
table1_df.to_csv(table1_handle, sep='\t', index=False)
table1_handle.close()