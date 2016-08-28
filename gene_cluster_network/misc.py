
import os
import re
import time

from Bio import Entrez

Entrez.email = "testing@ucsd.edu"


def parse_antiSMASH(content):
    """ Parse antiSMASH output
    """

    rule = r"""
        ^
        ClusterBlast\ scores\ for\ (?P<target>.*)\n+
        Table\ of\ genes,\ locations,\ strands\ and\ annotations\ of\ query\ cluster:\n+
        (?P<QueryCluster>
          (\w+ \s+ \d+ \s \d+ \s [+|-] \s \w+ (\s+\w+)?\s*\n+)+
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
            >>\n
            \d+ \. .+ \n
            Source:\ .+ \n
            Type:\ .+ \n
            Number\ of\ proteins\ with\ BLAST\ hits\ to\ this\ cluster:\ \d+ \n
            Cumulative\ BLAST\ score:\ \d+
            \n \n
            Table\ of\ genes,\ locations,\ strands\ and\ annotations\ of\ subject\ cluster:\n
            (?P<TableGenes>
              (\w+ \"? \t \w+ \t \d+ \t \d+ \t [+|-] \t .* \n)+
              )
            \n
            Table\ of\ Blast\ hits\ \(query\ gene,\ subject\ gene,\ %identity,\ blast\ score,\ %coverage,\ e-value\): \n
            (?P<BlastHit>
              (\w+ \t \w+ \"? \t \d+ \t \d+ \t \d+\.\d+ \t \d+\.\d+e[+|-]\d+ \t \n)+
              )
            \n+
            )+
          )
          \n*
          $
    """

    parsed = re.search(rule, content, re.VERBOSE).groupdict()

    output = {}

    output['target'] = parsed['target']

    output['QueryCluster'] = {'TableGenes': [], 'location_start': [],
            'location_end': [], 'strands': [], 
            'annotation': []}
    for row in re.finditer(r"""(\w+) \s+ (\d+) \s (\d+) \s ([+|-]) \s (\w+ (?:\s \w+)?) \s* \n+""", parsed['QueryCluster'], re.VERBOSE):
        output['QueryCluster']['TableGenes'].append(row.group(1))
        output['QueryCluster']['location_start'].append(int(row.group(2)))
        output['QueryCluster']['location_end'].append(int(row.group(3)))
        output['QueryCluster']['strands'].append(row.group(4))
        output['QueryCluster']['annotation'].append(row.group(5))

    #output['SignificantHits'] = {'id': [], 'name': [], 'description': []}
    output['SignificantHits'] = {}
    for row in re.finditer(r"""(?P<id>\d+) \. \ (?P<locus>\w+)_(?P<cluster>\w+) \t (?P<description>.*) \n+""", parsed['SignificantHits'], re.VERBOSE):
        hit = row.groupdict()

        if hit['locus'] not in output['SignificantHits']:
            output['SignificantHits'][hit['locus']] = {}

        if hit['cluster'] not in output['SignificantHits'][hit['locus']]:
            output['SignificantHits'][hit['locus']][hit['cluster']] = {}

        for v in ['id', 'description']:
            output['SignificantHits'][hit['locus']][hit['cluster']][v] = hit[v]

    for block in re.finditer(r"""
            >>\n
            (?P<id>\d+) \. \ (?P<locus>\w+)_(?P<cluster>\w+) \n
            Source:\ (?P<source>.+) \n
            Type:\ (?P<type>.+) \n
            Number\ of\ proteins\ with\ BLAST\ hits\ to\ this\ cluster:\ \d+ \n
            Cumulative\ BLAST\ score:\ \d+
            \n \n
            Table\ of\ genes,\ locations,\ strands\ and\ annotations\ of\ subject\ cluster:\n
            (?P<TableGenes>
              (\w+ \"? \t \w+ \t \d+ \t \d+ \t [+|-] \t .* \n)+
              )
            \n
            Table\ of\ Blast\ hits\ \(query\ gene,\ subject\ gene,\ %identity,\ blast\ score,\ %coverage,\ e-value\): \n
            (?P<BlastHit>
              (\w+ \t \w+ \"? \t \d+ \t \d+ \t \d+\.\d+ \t \d+\.\d+e[+|-]\d+ \t \n)+
              )
            \n+
""", parsed['Details'], re.VERBOSE):
        #block = block.groupdict()

        tmp = {'TableGenes': {}, 'TableBlast': {}}

        tmp = {
            'TableGenes': [], 'block': [], 'location_start': [],
            'location_end': [], 'strands': [], 'annotation': []}
        for row in re.finditer(r"""(\w+ \"?) \t (\w+) \t (\d+) \t (\d+) \t ([+|-]) \t (.*) \n""", block.groupdict()['TableGenes'], re.VERBOSE):
            tmp['TableGenes'].append(row.group(1))
            tmp['block'].append(row.group(2))
            tmp['location_start'].append(int(row.group(3)))
            tmp['location_end'].append(int(row.group(4)))
            tmp['strands'].append(row.group(5))
            tmp['annotation'].append(row.group(6))

        output['SignificantHits'][block.groupdict()['locus']][block.groupdict()['cluster']]['TableGenes'] = \
                tmp.copy()

        tmp = {
            'QueryGene': [], 'SubjectGene': [], 'Identity': [],
            'BlastScore': [], 'Coverage': [], 'e-value': []}
        for row in re.finditer(r"""(\w+) \t (\w+ \"?) \t (\d+) \t (\d+) \t (\d+(?:\.\d+)?) \t (\d+\.\d+e[+|-]\d+) \t \n""", block.groupdict()['BlastHit'], re.VERBOSE):
            tmp['QueryGene'].append(row.group(1))
            tmp['SubjectGene'].append(row.group(2))
            tmp['Identity'].append(row.group(3))
            tmp['BlastScore'].append(row.group(4))
            tmp['Coverage'].append(row.group(5))
            tmp['e-value'].append(row.group(6))

        output['SignificantHits'][block.groupdict()['locus']][block.groupdict()['cluster']]['TableBlast'] = \
                tmp.copy()

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

    for hit in c['SignificantHits'].keys():
        for cluster in c['SignificantHits'][hit]:
            table_genes = c['SignificantHits'][hit][cluster]['TableGenes']

            filename_out = os.path.join(
                    output_path,
                    "%s_%s_%s-%s.gbk" % (hit, cluster,
                        min(table_genes['location_start']),
                        max(table_genes['location_end'])))

            if os.path.isfile(filename_out):
                print "Already downloaded %s" % filename_out
            else:
                print "Requesting hit: %s, start: %s, end: %s" % (
                        hit,
                        min(table_genes['location_start']),
                        max(table_genes['location_end']))

                content = efetch_hit(
                    term=hit,
                    seq_start=min(table_genes['location_start']),
                    seq_stop=max(table_genes['location_end']))

                print "Saving %s" % filename_out
                with open(filename_out, 'w') as f:
                    f.write(content)
