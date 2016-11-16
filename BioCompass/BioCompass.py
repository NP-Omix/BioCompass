
import os
import re
import time
import json
from collections import OrderedDict
import pkg_resources

import pandas as pd
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "testing@ucsd.edu"


def untag(rule):
    return re.sub('\?P<.*?>', '', rule)

def parse_antiSMASH(content):
    """ Parse antiSMASH output
    """

    rule_table_genes = r"""
        (?P<subject_gene> \w+ \"?) \t
        \w+ \t
        (?P<location_start> \d+) \t
        (?P<location_end> \d+) \t
        (?P<strands> [+|-]) \t
        (?P<product> .*) \n
        """

    rule_table_blasthit = r"""
        (?P<query_gene> \w+ )\"? \t
        (?P<subject_gene> \w+ )\"? \t
        (?P<identity> \d+) \t
        (?P<blast_score> \d+) \t
        (?P<coverage> \d+(?:\.\d+)?) \t
        (?P<evalue> \d+\.\d+e[+|-]\d+) \t
        \n
        """

    rule_query_cluster = r"""
        (?P<query_gene> \w+) \s+
        (?P<location_start> \d+) \s
        (?P<location_end> \d+) \s
        (?P<strands> [+|-]) \s
        (?P<product> \w+ (?:\s \w+)?) \s* \n+
        """

    rule_detail = r"""
        >>\n
        (?P<id>\d+) \. \s+
            (?P<cluster_subject> (?P<locus>\w+)_(?P<cluster>\w+)) \n
        Source: \s+ (?P<source>.+?) \s* \n
        Type: \s+ (?P<type>.+) \s* \n
        Number\ of\ proteins\ with\ BLAST\ hits\ to\ this\ cluster:\ (?P<n_hits> \d+ ) \n
        Cumulative\ BLAST\ score:\ (?P<cum_BLAST_score> \d+ )
        \n \n
        Table\ of\ genes,\ locations,\ strands\ and\ annotations\ of\ subject\ cluster:\n
        (?P<TableGenes>
            (
            """ + untag(rule_table_genes) + r"""
            )+
        )
        \n
        Table\ of\ Blast\ hits\ \(query\ gene,\ subject\ gene,\ %identity,\ blast\ score,\ %coverage,\ e-value\): \n
        (?P<BlastHit>
            (\w+ \t \w+ \"? \t \d+ \t \d+ \t \d+\.\d+ \t \d+\.\d+e[+|-]\d+ \t \n)+
        )
        \n+
        """

    rule = r"""
        ^
        ClusterBlast\ scores\ for\ (?P<target>.*)\n+
        Table\ of\ genes,\ locations,\ strands\ and\ annotations\ of\ query\ cluster:\n+
        (?P<QueryCluster>
            (
            """ + untag(rule_query_cluster) + r"""
            )+
        )
        \n \n+
        Significant \  hits:\ \n
        (?P<SignificantHits>
          (\d+ \. \ \w+ \t .* \n+)+
          )
        \n \n
        (?P<Details>
          Details:\n\n
          (
            """ + untag(rule_detail) + r"""
          )+
        )
        \n*
        $
        """
    parsed = re.search(rule, content, re.VERBOSE).groupdict()

    output = {}
    for k in ['target', 'QueryCluster', 'SignificantHits']:
        output[k] = parsed[k]


    QueryCluster = OrderedDict()
    for k in re.search(
            rule_query_cluster, parsed['QueryCluster'],
            re.VERBOSE).groupdict().keys():
        QueryCluster[k] = []
    for row in re.finditer(
            rule_query_cluster, parsed['QueryCluster'], re.VERBOSE):
        row = row.groupdict()
        for k in row:
            QueryCluster[k].append(row[k])
    output['QueryCluster'] = QueryCluster


    output['SignificantHits'] = OrderedDict()
    for row in re.finditer(
            r"""(?P<id>\d+) \. \ (?P<cluster_subject> (?P<locus>\w+)_(?P<locus_cluster>\w+)) \t (?P<description>.*) \n+""", parsed['SignificantHits'], re.VERBOSE):
        hit = row.groupdict()
        cs = hit['cluster_subject']

        if cs not in output['SignificantHits']:
            output['SignificantHits'][cs] = OrderedDict()

        for v in ['id', 'description', 'locus', 'locus_cluster']:
            output['SignificantHits'][cs][v] = hit[v]

    for block in re.finditer(rule_detail, parsed['Details'], re.VERBOSE):
        block = dict(block.groupdict())

        content = block['TableGenes']
        block['TableGenes'] = OrderedDict()
        for k in re.findall('\(\?P<(.*?)>', rule_table_genes):
            block['TableGenes'][k] = []
        for row in re.finditer(rule_table_genes, content, re.VERBOSE):
            row = row.groupdict()
            for k in row:
                block['TableGenes'][k].append(row[k])

        content = block['BlastHit']
        block['BlastHit'] = OrderedDict()
        for k in re.findall('\(\?P<(.*?)>', rule_table_blasthit):
            block['BlastHit'][k] = []
        for row in re.finditer(rule_table_blasthit, content, re.VERBOSE):
            row = row.groupdict()
            for k in row:
                block['BlastHit'][k].append(row[k])

        for k in block:
            output['SignificantHits'][block['cluster_subject']][k] = block[k]

    return output


def antiSMASH_to_dataFrame(content):
    """ Extract an antiSMASH file as a pandas.DataFrame
    """
    parsed = parse_antiSMASH(content)
    output = pd.DataFrame()
    for cs in parsed['SignificantHits']:
        clusterSubject = parsed['SignificantHits'][cs].copy()
        df = pd.merge(
                pd.DataFrame(clusterSubject['BlastHit']),
                pd.DataFrame(clusterSubject['TableGenes']),
                on='subject_gene', how='outer')

        del(clusterSubject['BlastHit'])
        del(clusterSubject['TableGenes'])

        for v in clusterSubject:
            df[v] = clusterSubject[v]
        output = output.append(df, ignore_index=True)

    return output


class antiSMASH_file(object):
    """ A class to handle antiSMASH file output.
    """
    def __init__(self, filename):
        self.data = {}
        self.load(filename)

    def __getitem__(self, item):
        return self.data[item]

    def keys(self):
        return self.data.keys()

    def load(self, filename):
        self.data = {}
        with open(filename, 'r') as f:
            parsed = parse_antiSMASH(f.read())
            for v in parsed:
                self.data[v] = parsed[v]


def efetch_hit(term, seq_start, seq_stop):
    """ Fetch the relevant part of a hit
    """
    db = "nucleotide"

    maxtry = 3
    ntry = -1
    downloaded = False

    while ~downloaded and (ntry <= maxtry):
        ntry += 1
        try:
            handle = Entrez.esearch(db=db, term=term)
            record = Entrez.read(handle)

            assert len(record['IdList']) == 1, \
                    "Sorry, I'm not ready to handle more than one record"

            handle = Entrez.efetch(db=db, rettype="gb", retmode="text",
                    id=record['IdList'][0],
                    seq_start=seq_start, seq_stop=seq_stop)
            content = handle.read()
            downloaded = True
        except:
            nap = ntry*3
            print "Fail to download (term). I'll take a nap of %s seconds ", \
                    " and try again."
            time.sleep(ntry*3)

    return content


def download_hits(filename, output_path):
    """ Download the GenBank block for all hits by antiSMASH
    """
    c = antiSMASH_file(filename)

    for cs in c['SignificantHits'].keys():
        locus = c['SignificantHits'][cs]['locus']
        table_genes = c['SignificantHits'][cs]['TableGenes']

        filename_out = os.path.join(
                output_path,
                "%s_%s-%s.gbk" % (locus,
                    min(table_genes['location_start']),
                    max(table_genes['location_end'])))

        if os.path.isfile(filename_out):
            print "Already downloaded %s" % filename_out
        else:
            print "Requesting cluster_subject: %s, start: %s, end: %s" % (
                    locus,
                    min(table_genes['location_start']),
                    max(table_genes['location_end']))

            content = efetch_hit(
                term=locus,
                seq_start=min(table_genes['location_start']),
                seq_stop=max(table_genes['location_end']))

            print "Saving %s" % filename_out
            with open(filename_out, 'w') as f:
                f.write(content)


import urlparse
import urllib2
import tempfile
import tarfile
import os
def download_mibig(outputdir, version='1.3'):
    """ Download and extract MIBiG files into outputdir
    """
    assert version in ['1.0', '1.1', '1.2', '1.3'], \
            "Invalid version of MIBiG"

    server = 'http://mibig.secondarymetabolites.org'
    filename = "mibig_gbk_%s.tar.gz" % version
    url = urlparse.urljoin(server, filename)

    with tempfile.NamedTemporaryFile(delete=True) as f:
        u = urllib2.urlopen(url)
        f.write(u.read())
        f.file.flush()
        tar = tarfile.open(f.name)
        tar.extractall(path=outputdir)
        tar.close()

    # MIBiG was packed with strange files ._*gbk. Let's remove it
    for f in [f for f in os.listdir(outputdir) if f[:2] == '._']:
        os.remove(os.path.join(outputdir, f))

#def gbk2tablegen(gb_file, strain_id=None):
#def cds_from_gbk(gb_file, strain_id=None):
def cds_from_gbk(gb_file):
    gb_record = SeqIO.read(open(gb_file,"rU"), "genbank")

    #if strain_id is not None:
    #    gb_record.id = strain_id

    output = pd.DataFrame()
    sign = lambda x: '+' if x > 0 else '-'
    for feature in gb_record.features:
        if feature.type == "CDS":
            tmp = {}
            tmp = {'BGC': gb_record.id,
                    'locus_tag': feature.qualifiers['locus_tag'][0],
                    'start': feature.location.start.position,
                    'stop': feature.location.end.position,
                    'strand': sign(feature.location.strand) }
            if 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    product = re.search( r"""smCOG: \s (?P<product>.*?) \s+ \(Score: \s* (?P<score>.*); \s* E-value: \s (?P<e_value>.*?)\);""", note, re.VERBOSE)
                    if product is not None:
                        product = product.groupdict()
                        product['score'] = float(product['score'])
                        product['e_value'] = float(product['e_value'])
                        for p in product:
                            tmp[p] = product[p]
            output = output.append(pd.Series(tmp), ignore_index=True)
    return output


def find_category_from_product(df):
    subcluster = json.loads(
            pkg_resources.resource_string(
                __name__, 'subcluster_dictionary.json'))
    def get_category(product):
        for s in subcluster:
            if re.search(s, product):
                return subcluster[s]
        return 'hypothetical'

    idx = df['product'].notnull()
    df['category'] = df.loc[idx, 'product'].apply(get_category)
    df['category'].fillna('hypothetical', inplace=True)
    return df


def get_hits(filename, criteria='cum_BLAST_score'):
    """

        Reproduces original Tiago's code: table_1_extender.py

        In the future allow different criteria. Right now it takes
          from the very first block, which has the highest Cumulative
          BLAST.
    """
    with open(filename) as f:
        df = antiSMASH_to_dataFrame(f.read())

    df.dropna(subset=['query_gene'], inplace=True)
    df.sort_values(by=criteria, ascending=False, na_position='last',
            inplace=True)
    return df.groupby('query_gene', as_index=False).first()
