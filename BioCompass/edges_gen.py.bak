import pandas as pd
import re
import os.path
import itertools
from sys import argv

script, cluster_name, itineration = argv

table2_df = pd.read_csv('../tables/%s_table2_%s.csv'%(cluster_name,itineration), sep='\t')

col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
col6 = []
col7 = []
col8 = []
col9 = []

for index, row in table2_df.iterrows():
    file_name = './%s_%s_%s/clusterblast_output.txt'%(table2_df.BGC.loc[index],itineration,table2_df.subcluster.loc[index])
    if os.path.isfile(file_name):
        with open(file_name).xreadlines() as f:
            for line in f:
                hit = re.search(r'^(\d*). (\S*)_(\S*)$',line)
                matching_genes = re.search(r'^Number of proteins with BLAST hits to this cluster: (\d*)$', line)
                mgb_score = re.search(r'^MultiGeneBlast score: (\d*).(\d*)$', line)
                blast_score = re.search(r'^Cumulative Blast bit score: (\d*)$',line)
                if hit != None:
                    col1.append('%s'%table2_df.BGC.loc[index])
                    col2.append('%s'%table2_df.subcluster.loc[index])
                    col3.append('%s'%table2_df.category.loc[index])
                    col4.append('%s'%table2_df.CDSs.loc[index])
                    col5.append(hit.group(2))
                    col9.append('%s'%itineration)
                if mgb_score != None:
                    col6.append('%s.%s'%mgb_score.group(1,2))
                if blast_score != None:
                    col7.append(blast_score.group(1))
                if matching_genes != None:
                    perc_matching = int(matching_genes.group(1))*100/(table2_df.CDSs.loc[index])
                    col8.append(perc_matching)
                
    else:
        col1.append('%s'%table2_df.BGC.loc[index])
        col2.append('%s'%table2_df.subcluster.loc[index])
        col3.append('%s'%table2_df.category.loc[index])
        col4.append('%s'%table2_df.CDSs.loc[index])
        col5.append('%s'%table2_df.BGC.loc[index])
        col6.append('NA')
        col7.append('NA')
        col8.append('NA')
        col9.append('%s'%itineration)
        
frames = {'BGC':col1, 'BGC_subcluster':col2, 'BGC_subcluster_category':col3, 'BGC_subcluster_genes':col4, 'BLAST_hit':col5, 'MultiGeneBlast_score':col6, 'BLAST_score':col7, 'percent_matching_genes':col8, 'BGC_iteration':col9}

edges_df = pd.DataFrame(frames)

edges_handle = open('%s_edges_%s.txt'%(cluster_name,itineration), "w")
edges_df.to_csv(edges_handle, sep='\t', index=False)
edges_handle.close()

