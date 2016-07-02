from Bio import SeqIO
import pandas as pd
import re
from Bio.SeqUtils import GC
from sys import argv

script, strain_name, edges_file = argv

edges_df = pd.read_csv(edges_file, sep='\t')

ref_list = edges_df['BGC'].drop_duplicates(inplace=False)

hits_list = edges_df['BLAST_hit'].drop_duplicates(inplace=False)


col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
col6 = []
col7 = []
col8 = []


for item in ref_list:
    gb_record = SeqIO.read(open('../tables/%s.gbk'%item,"r"), "genbank")
    table_1 = pd.read_csv('../tables/%s_table1.csv'%item, sep='\t')
    for feature in gb_record.features:
            if feature.type == "cluster":
                col1.append(item)
                col2.append(feature.qualifiers['product'])
                col3.append('%s'%strain_name)
                col4.append(gb_record.description)
                col5.append(len(table_1))
                size = re.search(r'(\d*):(\d*)(.*)','%s'%feature.location)
                col6.append(size.group(2))
                GC_cont = GC(str(gb_record.seq))
                col7.append(round(GC_cont,2))
                if 'N' in gb_record.seq:
                    col8.append('incomplete')
                else:
                    col8.append('complete')
                        
for item in hits_list:
    l = re.search(r'^%s(.*)$'%strain_name, str(item))
    m = re.search(r'^BGC(.*)$', str(item))
    n = re.search(r'(\S*)_\d\d\d', str(item))
    if n != None and l == None:
        gb_record = SeqIO.read(open('../../database_clusters/%s/%s.gbk'%(n.group(1),item),"r"), "genbank")
        for feature in gb_record.features:
            if feature.type == "cluster":
                col1.append(item)
                col2.append(feature.qualifiers['product'])
                col3.append('%s'%n.group(1))
                col4.append(gb_record.description)
                CDSs = []
                for feature in gb_record.features:
                    if feature.type == 'CDS':
                        CDSs.append(feature)
                col5.append(len(CDSs))
                col6.append(len(gb_record.seq))
                GC_cont = GC(str(gb_record.seq))
                col7.append(round(GC_cont,2))
                if 'N' in gb_record.seq:
                    col8.append('incomplete')
                else:
                    col8.append('complete')
    if m != None:
        gb_record = SeqIO.read(open('../../database_clusters/MIBiG/%s.gbk'%item,"r"), "genbank")
        col1.append(item)
        col2.append('MIBiG_%s'%gb_record.description)
        col3.append('MIBiG')
        col4.append(gb_record.annotations['organism'])
        CDSs = []
        for feature in gb_record.features:
            if feature.type == 'CDS':
                CDSs.append(feature)
        col5.append(len(CDSs))
        col6.append(len(gb_record.seq))
        GC_cont = GC(str(gb_record.seq))
        col7.append(round(GC_cont,2))
        if 'N' in gb_record.seq:
            col8.append('incomplete')
        else:
            col8.append('complete')
        
frames = {'BGC':col1,'Category':col2, 'Origin':col3, 'Organism':col4, 'Number_of_genes':col5, 'Size(bp)':col6, 'GC_content':col7, 'Completness':col8}

features_df = pd.DataFrame(frames, index=None)

features_handle = open('%s_features.txt' % strain_name, "w")
features_df.to_csv(features_handle, sep='\t', index=False)
features_handle.close()

