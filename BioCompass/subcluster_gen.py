from collections import defaultdict
from sys import argv
import pandas as pd
import re
import string
import numpy as np
from sklearn.cluster import DBSCAN
from collections import Counter

script, strain_name = argv

table1_df = pd.read_csv('%s_table1.csv' % strain_name, sep='\t')
table1_df['product'].fillna('None', inplace=True)

#This first portion will create the distance matrix

def score_match(table, index_gene1, index_gene2):
    score = 0
    gene1_category = table1_df.category.loc[index_gene1]
    gene2_category = table1_df.category.loc[index_gene2]
    if  gene1_category == 'hypothetical' or gene2_category == 'hypothetical':
        if gene1_category == gene2_category:
            score = score + 1
        elif gene1_category != gene2_category:
            score = score + 2
    else:
        if gene1_category == gene2_category:
            score = score + 5
    gene1_best_BGC = table1_df.best_hit_BGC.loc[index_gene1]
    gene2_best_BGC = table1_df.best_hit_BGC.loc[index_gene2]
    if gene1_best_BGC != 'None' and gene2_best_BGC != 'None':
        if gene1_best_BGC == gene2_best_BGC:
            score = score + 2
            gene1_best_hit_pos = re.search(r'^\D*([0-9]*)',table1_df.best_hit_gene_loc.loc[index_gene1])
            gene2_best_hit_pos = re.search(r'^\D*([0-9]*)',table1_df.best_hit_gene_loc.loc[index_gene2])
            dif_best_hit_pos = abs(abs((int(gene2_best_hit_pos.group(1)) - int(gene1_best_hit_pos.group(1)))) - abs((index_gene2-index_gene1)))
            if dif_best_hit_pos == 0:
                score = score + 3
            elif dif_best_hit_pos == 1:
                score = score + 2
            elif dif_best_hit_pos == 2:
                score = score + 1
    else:
        score = score + 1
    return score

for index,row in table1_df.iterrows():
    scores = []
    for gene in range(0,len(table1_df)):
        scores.append(score_match(table1_df,gene,index))
    if index == 0:
        A = np.vstack([scores])
    else:
        A = np.vstack([A,scores])

#This second portion will run dbscan to create a subclusters possibilities

def repeated(db_arrays,db):
    for i in range(0,len(db_arrays)-1):
        if np.array_equal(db_arrays[i],db) == False:
            continue
        else:
            return True
            break

def parse_db(db):
    D = defaultdict(list)
    for i,item in enumerate(db):
        D[item].append(i)
    D = {k:v for k,v in D.items() if len(v)>1}
    return D

def find_category(categories,col5):
    if 'biosynthetic' in categories:
        col5.append('biosynthetic')
    else:
        if len(categories) > 1:
            category = re.search(r'^\[\(\'(\S*)\'',str(Counter(categories).most_common(1))).group(1)
            col5.append('%s'%category)
        else:
            col5.append('%s'%categories[0])
    
count = 0
    
for itn in range(1,len(A)):
    db = DBSCAN(eps=itn, min_samples=2).fit_predict(A)
    if itn == 1:
        db_arrays = np.vstack([db])
    else:
        db_arrays = np.vstack([db_arrays,db])
    if repeated(db_arrays,db) == True:
        continue
    else:
        subcluster_dict = parse_db(db)
        col1 = []
        col2 = []
        col3 = []
        col4 = []
        col5 = []
        for key, value in subcluster_dict.iteritems():
            col1.append(strain_name)
            col2.append(string.ascii_uppercase[subcluster_dict.keys().index(key)])
            col3.append(len(value))
            categories = []
            genes = []
            for item in value:
                categories.append(table1_df.category.loc[item])
                genes.append(table1_df.locus_tag.loc[item])
            genes = ','.join(genes)
            col4.append(genes)
            find_category(categories,col5)
        frames = {'BGC':col1, 'subcluster':col2, 'CDSs':col3, 'loci':col4, 'category':col5}
        count = count + 1
        table2_df = pd.DataFrame(frames, index=None)
        table2_df.to_csv('%s_table2_%d.csv' % (strain_name,count), sep='\t', index=False)
