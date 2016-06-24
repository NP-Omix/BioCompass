from __future__ import division
from sys import argv
import pandas as pd
import re
import string
import collections
from collections import Counter
import numpy as np
from sklearn.cluster import DBSCAN

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
            gene1_best_hit_pos = re.search(r'^(\S*)_(\d*)',table1_df.best_hit_gene.loc[index_gene1])
            gene2_best_hit_pos = re.search(r'^(\S*)_(\d*)',table1_df.best_hit_gene.loc[index_gene2])
            dif_best_hit_pos = abs(abs((int(gene2_best_hit_pos.group(2)) - int(gene1_best_hit_pos.group(2)))) - abs((index_gene2-index_gene1)))
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

def starter_line(i):
    col1.append(strain_name)
    col2.append(string.ascii_uppercase[0])
    col3.append(1)
    col4.append(table1_df.start.loc[i])
    
def continuing(i):
    col3[len(col3)-1] = col3[len(col3)-1] + 1
    categories.append(table1_df.category.loc[i])
    
    
def find_category(categories,colpre6):
    if 'biosynthetic' in categories:
        colpre6.append('biosynthetic')
    else:
        if len(categories) > 1:
            category = re.search(r'^\[\(\'(\S*)\'',str(Counter(categories).most_common(1))).group(1)
            colpre6.append('%s'%category)
        else:
            colpre6.append('%s'%categories[0])

def new_subcluster(i):
    col1.append(strain_name)
    if len(col2) < 26:
        col2.append(string.ascii_uppercase[len(col2)])
    else:
        col2.append('A%s'%string.ascii_uppercase[len(col2)-26])
    col3.append(1)
    col4.append(table1_df.start.loc[i])
    col5.append(table1_df.stop.loc[i-1])
    find_category(categories,colpre6) 

def col6_gen(colpre6):
    seen = collections.defaultdict(list)
    for item in colpre6:
        if item not in seen:
            seen[item] = 1
            col6.append(item+"_%d"%seen[item])
        else:
            m = seen[item]
            seen[item] = m + 1
            col6.append(item+"_%d"%seen[item])

def repeated(db_arrays,db):
    for i in range(0,len(db_arrays)-1):
        if np.array_equal(db_arrays[i],db) == False:
            continue
        else:
            return True
            break

count = 0

for index in range(1,len(A)):
    db = DBSCAN(eps=index, min_samples=2).fit_predict(A)
    if index == 1:
        db_arrays = np.vstack([db])
    else:
        db_arrays = np.vstack([db_arrays,db])
    if repeated(db_arrays,db) == True:
        continue
    else:
        col1 = []
        col2 = []
        col3 = []
        col4 = []
        col5 = []
        colpre6 = []
        col6 = []
        for i in range(0,len(table1_df)):
            if i == 0:
                starter_line(1)
                categories = []
                categories.append(table1_df.category.loc[i])
            if i > 0 and i < len(table1_df)-1:
                if db[i] == db[i-1]:
                    continuing(i)
                else:
                    if db[i-1] == db[i+1] and table1_df.category.loc[i] != 'biosynthetic':
                        continuing(i)
                    elif db[i-2] == db[i]:
                        continuing(i)
                    else:
                        new_subcluster(i)
                        categories = []
                        categories.append(table1_df.category.loc[i])
            if i == len(table1_df)-1:
                if db[i] == db[i-1]:
                    col3[len(col3)-1] = col3[len(col3)-1] + 1
                    col5.append(table1_df.stop.loc[i])
                    categories.append(table1_df.category.loc[i])
                    find_category(categories,colpre6)
                else:
                    if db[i-2] == db[i]:
                        col3[len(col3)-1] = col3[len(col3)-1] + 1
                        col5.append(table1_df.stop.loc[i-1]) 
                        categories.append(table1_df.category.loc[i])
                        find_category(categories,colpre6)
                    elif table1_df.category.loc[i] == 'hypothetical':
                        col3[len(col3)-1] = col3[len(col3)-1] + 1
                        col5.append(table1_df.stop.loc[i-1])
                        find_category(categories,colpre6)
                    else:
                        new_subcluster(i)
                        col5.append(table1_df.stop.loc[i])
                        colpre6.append(table1_df.category.loc[i])
        col6_gen(colpre6)
        frames = {'BGC':col1, 'subcluster':col2, 'CDSs':col3, 'start':col4, 'stop':col5, 'category':col6}
        count = count + 1
        table2_df = pd.DataFrame(frames, index=None)
        table2_df.to_csv('%s_table2_%d.csv' % (strain_name,count), sep='\t', index=False)
