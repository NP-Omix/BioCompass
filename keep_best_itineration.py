import pandas as pd
import itertools
from sys import argv

script, strain_name = argv

edges_df = pd.read_csv('%s_edges_all_itineration.txt'%strain_name, sep='\t')

edges_df.fillna(0, inplace=True)

edges_df_filtered = edges_df.drop_duplicates(subset=['BGC','BGC_iteration','BGC_subcluster'],inplace=False)

edges_df_filtered.reset_index(inplace=True,drop=True)

col1 = []
col2 = []
col3 = []
col4 = []

for index,row in edges_df_filtered.iterrows():
    if index == 0:
        col1.append(edges_df_filtered.BGC.loc[index])
        col2.append(edges_df_filtered.BGC_iteration.loc[index])
        col3.append(1)
        col4.append(edges_df_filtered.BLAST_score.loc[index])
    if index > 0 and index < len(edges_df_filtered):
        if edges_df_filtered.BGC.loc[index] == edges_df_filtered.BGC.loc[index-1] and edges_df_filtered.BGC_iteration.loc[index] == edges_df_filtered.BGC_iteration.loc[index-1]:
            col3[len(col3)-1] = col3[len(col3)-1] + 1
            col4[len(col3)-1] = col4[len(col4)-1] + edges_df_filtered.BLAST_score.loc[index]
        else:
            col1.append(edges_df_filtered.BGC.loc[index])
            col2.append(edges_df_filtered.BGC_iteration.loc[index])
            col3.append(1)
            col4.append(edges_df_filtered.BLAST_score.loc[index])

frames = {'BGC':col1, 'Itineration':col2, 'Number_of_subclusters':col3, 'Total_score':col4}

table3_df = pd.DataFrame(frames)

table3_df.reset_index(inplace=True,drop=True)

table3_handle = open('%s_table3.csv' % strain_name, "w")
table3_df.to_csv(table3_handle, sep='\t', index=False)
table3_handle.close()

table3_df_sorted = table3_df.sort_values(by=['BGC','Number_of_subclusters'],ascending=False,inplace=False)

table3_df_sorted.sort_values(by=['BGC','Total_score'],ascending=True,inplace=True)

table3_df_sorted.drop_duplicates(subset=['BGC'],inplace=True,keep='last')

table3_dict = {}
indexes_to_drop = []

for i,r in table3_df_sorted.iterrows():
    table3_dict[table3_df_sorted.BGC.loc[i]] = table3_df_sorted.Itineration.loc[i]
    
for index,row in edges_df.iterrows():
    if edges_df.BGC.loc[index] in table3_dict:
        if edges_df.BGC_iteration.loc[index] != table3_dict[edges_df.BGC.loc[index]]:
            indexes_to_drop.append(index)

edges_final = edges_df.drop(edges_df.index[indexes_to_drop])

edges_final.replace(to_replace=0, inplace=True, value='NA')

edges_final_handle = open('%s_edges_best_itineration.txt' % strain_name, "w")
edges_final.to_csv(edges_final_handle, sep='\t', index=False)
edges_final_handle.close()