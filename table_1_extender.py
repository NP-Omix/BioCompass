from sys import argv
from Bio import SeqIO
import pandas as pd
import re
import itertools
import os.path

script, strain_name = argv

table1_df = []

table1_df = pd.read_csv('../outputs/tables/%s_table1.csv'%strain_name, sep='\t')

def find_boarders(file_name):
    with open(file_name).xreadlines() as f:
        for num, line in enumerate(f):
            header = re.search(r'^(\d*). (\S*)_\d$',line)
            if header != None and num not in numbers:
                numbers.append(num)
                headers.append(header.group())
    temp_dict = dict(zip(numbers, headers))
    return temp_dict
                
def itinerate_temp_file(file_name,temp_dict):
    with open(file_name) as fp:
        temp_file = open("temp.txt", "w")
        for num, line in enumerate(fp):
            if i[0] != sorted(temp_dict)[-1]:
                if num >= i[0] and num < j[0]:
                    temp_file.writelines(line)
            else:
                if num >= i[0]:
                    temp_file.writelines(line)
        temp_file.close()

    
def find_best_hits(table1_df):
    for item in table1_df['locus_tag']:
        with open("temp.txt", "r").xreadlines() as tf:
            for line in tf:
                best_hit = re.search(r'^%s\t(\S*)\t(\d*)\t(\d*)\t(\d*)'%item,line)
                if best_hit and item not in col7:
                    col7.append(item)
                    hit = re.search(r'^(\d*). (\S*)_\d$',i[1])
                    col8.append(hit.group(2))
                    col9.append(best_hit.group(1))
                    col10.append(best_hit.group(2))
                    col11.append(best_hit.group(4))

#find boarders for temp_file
file_name = '../outputs/pre_mgb_result/%s/clusterblast_output.txt' % strain_name
numbers = []
headers = []
temp_dict = []
if os.path.isfile(file_name):
    temp_dict = find_boarders(file_name)
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []
    it = iter(sorted(temp_dict.iteritems()))
    i = it.next()
    if len(temp_dict) > 1:
        j = it.next()
    while True:
        if i[0] == sorted(temp_dict)[-1]:
            itinerate_temp_file(file_name,temp_dict)
            find_best_hits(table1_df)
            break
        else:
            itinerate_temp_file(file_name,temp_dict)
            find_best_hits(table1_df)
            i = j
            if j[0] != sorted(temp_dict)[-1]:
                j = it.next()
    frames = {'locus_tag':col7,'best_hit_BGC':col8,'best_hit_gene':col9,'best_hit_%id':col10,'best_hit_%cov':col11}
    new_cols_df = pd.DataFrame(frames, index=None)
    table1_df = pd.merge(table1_df, new_cols_df, on='locus_tag', how='outer')
    table1_df.fillna('None', inplace=True)       
else:
    col7 = []
    col8 = []
    for item in table1_df['locus_tag']:
        col7.append(item)
        col8.append('None')
    frames = {'locus_tag':col7,'best_hit_BGC':col8}
    new_cols_df = pd.DataFrame(frames, index=None)
    table1_df = pd.merge(table1_df, new_cols_df, on='locus_tag', how='outer')

table1_handle = open('../outputs/tables/%s_table1.csv' % strain_name, "w")
table1_df.to_csv(table1_handle, sep='\t', index=False)
table1_handle.close()