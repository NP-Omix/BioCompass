from sys import argv
from Bio import SeqIO
import pandas as pd
import re
import itertools
import os.path

script, strain_name, cluster_number = argv

table1_df = pd.read_csv('../outputs/tables/%s_%s_table1.csv'%(strain_name,cluster_number), sep='\t')

def find_boarders(file_name):
    with open(file_name).xreadlines() as f:
        for num, line in enumerate(f):
            header = re.search(r'^(\d*). (\S*)_(\S*)$',line)
            if header != None and num not in numbers:
                numbers.append(num)
                headers.append(header.group())
    temp_dict = dict(zip(numbers, headers))
    return temp_dict
                
def itinerate_temp_file(file_name,temp_dict):
    """ Extract from file_name a section defined by i into temp.txt"""
    with open(file_name) as fp:
        temp_file = open("temp.txt", "w")
        for num, line in enumerate(fp):
            # Is not the last one
            if i[0] != sorted(temp_dict)[-1]:
                if num >= i[0] and num < j[0]:
                    temp_file.writelines(line)
            else:
                if num >= i[0]:
                    temp_file.writelines(line)
        temp_file.close()

    
def find_best_hits(table1_df):
    """ For each ctg1_* @ table1_df['locus_tag'], get
    """
    for item in table1_df['locus_tag']:
        with open("temp.txt", "r").xreadlines() as tf:
            for line in tf:
                # query gene, subject gene, %identity, blast score, %coverage, e-value
                # e.g line <- 'ctg1_160\tSCATT_35940\t69\t169\t100.0\t2.3e-40\t\n'
                best_hit = re.search(r'^%s\t(\S*)\t(\d*)\t(\d*)\t(\d*)'%item,line)
                if best_hit and item not in col7:
                    # e.g. ctg1_160
                    col7.append(item)
                    # e.g i[1] <- '1. CP003219_c13'
                    hit = re.search(r'^(\d*). (\S*)_(\S*)',i[1])
                    # e.g. 'CP003219'
                    col8.append(hit.group(2))
                    # e.g. 'SCATT_35940'
                    col9.append(best_hit.group(1))
                    # e.g. '69'
                    col10.append(best_hit.group(2))
                    # e.g. '100'
                    col11.append(best_hit.group(4))

#find boarders for temp_file
short_cluster_number = re.search(r'0*([0-9]*)',cluster_number).group(1)
file_name = '../antiSMASH_input/%s/clusterblast/cluster%s.txt' % (strain_name,short_cluster_number)
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
    # e.g. i <- (139, '1. CP003219_c13')
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
    """
    ctg1_160        SCATT_35940     69      169     100.0   2.3e-40
    hence:
    best_hit_loc ->
      SCATT_35940     AEW95965        3878018 3878386 +       50S_ribosomal_protein_L14
    """
    col12 = []
    for item in col9:
        with open(file_name, "r").xreadlines() as tf:
            seen = []
            for line in tf:
                best_hit_loc = re.search(r'^%s\t(\S*)\t(.*)'%item,line)
                if best_hit_loc and item not in seen:
                    col12.append(best_hit_loc.group(1))
                    seen.append(item)
    frames = {'locus_tag':col7,'best_hit_BGC':col8,'best_hit_gene':col9,'best_hit_%id':col10,'best_hit_%cov':col11,'best_hit_gene_loc':col12}
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


table1_handle = pd.HDFStore('../outputs/tables/table1.h5')
table1_handle['%s_%s' % (strain_name,cluster_number)] = table1_df
table1_handle.close()


table1_handle = open('../outputs/tables/%s_%s_table1.csv'%(strain_name,cluster_number), "w")
table1_df.to_csv(table1_handle, sep='\t', index=False)
table1_handle.close()
