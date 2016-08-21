
import re


def parse_antiSMASH(content):
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

    output['SignificantHits'] = {'id': [], 'name': [], 'description': []}
    for row in re.finditer(r"""(\d+) \. \ (\w+) \t (.*) \n+""", parsed['SignificantHits'], re.VERBOSE):
        output['SignificantHits']['id'].append(row.group(1))
        output['SignificantHits']['name'].append(row.group(2))
        output['SignificantHits']['description'].append(row.group(3))

    output['Details'] = {}

    for block in re.finditer(r"""
            >>\n
            (?P<id>\d+) \. \ (?P<name>.+) \n
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
        name = block.groupdict()['name']
        output['Details'][name] = {'TableGenes': {}, 'TableBlast': {}}

        output['Details'][name]['TableGenes'] = {
            'TableGenes': [], 'block': [], 'location_start': [],
            'location_end': [], 'strands': [], 'annotation': []}
        for row in re.finditer(r"""(\w+ \"?) \t (\w+) \t (\d+) \t (\d+) \t ([+|-]) \t (.*) \n""", block.groupdict()['TableGenes'], re.VERBOSE):
            output['Details'][name]['TableGenes']['TableGenes'].append(row.group(1))
            output['Details'][name]['TableGenes']['block'].append(row.group(2))
            output['Details'][name]['TableGenes']['location_start'].append(row.group(3))
            output['Details'][name]['TableGenes']['location_end'].append(row.group(4))
            output['Details'][name]['TableGenes']['strands'].append(row.group(5))
            output['Details'][name]['TableGenes']['annotation'].append(row.group(6))

        output['Details'][name]['TableBlast'] = {
            'QueryGene': [], 'SubjectGene': [], 'Identity': [],
            'BlastScore': [], 'Coverage': [], 'e-value': []}
        for row in re.finditer(r"""(\w+) \t (\w+ \"?) \t (\d+) \t (\d+) \t (\d+(?:\.\d+)?) \t (\d+\.\d+e[+|-]\d+) \t \n""", block.groupdict()['BlastHit'], re.VERBOSE):
            output['Details'][name]['TableBlast']['QueryGene'].append(row.group(1))
            output['Details'][name]['TableBlast']['SubjectGene'].append(row.group(2))
            output['Details'][name]['TableBlast']['Identity'].append(row.group(3))
            output['Details'][name]['TableBlast']['BlastScore'].append(row.group(4))
            output['Details'][name]['TableBlast']['Coverage'].append(row.group(5))
            output['Details'][name]['TableBlast']['e-value'].append(row.group(6))

    return output
