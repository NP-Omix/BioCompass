import os.path
from sys import argv
from Bio import SeqIO
import pandas as pd
import re

script, gb_file, strain_name = argv

strain_id = os.path.basename(strain_name)

gb_record = SeqIO.read(open(gb_file,"r"), "genbank")

col2 = []
col3 = []
col4 = []
col5 = []
col6 = []

for (index, feature) in enumerate(gb_record.features):
    if feature.type == "CDS":
        locus_tag = re.search(r'^\[\'(\S*)\'\]$', '%s'%feature.qualifiers['locus_tag'])
        col2.append(locus_tag.group(1))
        col3.append(feature.location.start)
        col4.append(feature.location.end)
        if 'note' in feature.qualifiers:
            smCOG = re.search(r'smCOG:\s(.*)\s\(Score: (.*)', '%s'%feature.qualifiers['note'], re.S)
            if smCOG:
                col5.append(smCOG.group(1))
            else:
                col5.append(None)
        else:
            col5.append(None)
        col6.append(feature.location.strand)

col1 = [strain_id] * len(col2)
frames = {'BGC':col1, 'locus_tag':col2, 'start':col3, 'stop':col4, 'product':col5, 'strand':col6}

table1_df = pd.DataFrame(frames, index=None)

table1_handle = open('%s_table1.csv' % strain_name, "w")
table1_df.to_csv(table1_handle, sep='\t', index=False)
table1_handle.close()

gb_record.id = '%s' % strain_id

output_handle = open("%s.gbk"%strain_name, "w")
SeqIO.write(gb_record, output_handle, "genbank")
output_handle.close()
