import pandas as pd
import itertools
from sys import argv

script, strain_name = argv

edges_final = pd.read_csv('%s_edges_best_itineration.txt'%strain_name, sep='\t')


while True:
    try:
        per_genes = int(raw_input("\n What minimum percentage of matching genes would you like [from 0 to 100]? "))
        blast_score = int(raw_input("\n What minimum cumulative BLAST score would you like [from %s to %s]? "%(edges_final.BLAST_score.min(),edges_final.BLAST_score.max())))
        mgb_score = int(raw_input("\n What minimum multigeneBLAST score [from %s to %s]? "%(edges_final.MultiGeneBlast_score.min(),edges_final.MultiGeneBlast_score.max())))
        break
    except ValueError:
        print("Sorry, I didn't understand that. Try again using INTERGERS ONLY.")
        continue
    else:
        break
        
seen = []
indexes_to_drop = []
        
for index,row in edges_final.iterrows():
    query = '%s_%s'%(edges_final.BGC.loc[index],edges_final.BGC_subcluster.loc[index])
    if edges_final.percent_matching_genes.loc[index] >= per_genes and edges_final.BLAST_score.loc[index] >= blast_score and edges_final.MultiGeneBlast_score.loc[index] >= mgb_score:
        seen.append(query)
        continue
    else:
        if query in seen:
            indexes_to_drop.append(index)
        else:
            seen.append(query)
            edges_final.set_value(index,'BLAST_hit','%s'%edges_final.BGC.loc[index])
            edges_final.set_value(index,'percent_matching_genes',None)
            edges_final.set_value(index,'BLAST_score',None)
            edges_final.set_value(index,'MultiGeneBlast_score',None)

edges_final = edges_final.drop(edges_final.index[indexes_to_drop])
            
while True:
    self_loop = raw_input("\n Would you like to include subclusters with no matches?[Yes or No only!] ")
    if self_loop == "Yes" or self_loop == "yes":
        edges_final.MultiGeneBlast_score.replace(to_replace='NaN', value=mgb_score, inplace=True)
        edges_final.BLAST_score.replace(to_replace='NaN', value=blast_score, inplace=True)
        edges_final.percent_matching_genes.replace(to_replace='NaN', value=100, inplace=True)
        break
    elif self_loop == "No" or self_loop == "no":
        break
    else:
        print("Sorry, answer Yes or No only")
        continue

edges_final_handle = open('%s_edges_filtered.txt'%strain_name, "w")
edges_final.to_csv(edges_final_handle, sep='\t', index=False)
edges_final_handle.close()