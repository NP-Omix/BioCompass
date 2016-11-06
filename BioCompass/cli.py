# -*- coding: utf-8 -*-

import click

from .BioCompass import download_hits
from .BioCompass import download_mibig

@click.group()
def main():
    pass

@main.command(name="download-hits")
@click.option('--outputdir', default='./', type=click.Path(exists=True),
        help="Path to save the NCBI clusters.")
@click.argument('mgbfile', type=click.Path(exists=True))
        #help="Multigeneblast file containing NCBI references to be downloaded.")
def downloadHits(mgbfile, outputdir):
    """Download NCBI clusters listed t in multigeneblast file."""
    download_hits(mgbfile, outputdir)

@main.command(name="download-MIBiG")
@click.option('--outputdir', default='./', type=click.Path(exists=True),
        help="Path to save the MIBig genbank files.")
@click.option('--version', type=unicode, default='1.3',
        help="Version of MIBiG to download.")
def downloadMIBiG(outputdir, version):
    """Download MIBiG gbk database."""
    download_mibig(outputdir, version=version)

if __name__ == "__main__":
    main()
