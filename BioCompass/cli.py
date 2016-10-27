# -*- coding: utf-8 -*-

import click

from BioCompass.misc import download_hits

@click.command()
@click.option('--outputdir', default='./', type=click.Path(exists=True),
        help="Path to save the NCBI clusters.")
@click.argument('mgbfile', type=click.Path(exists=True))
        #help="Multigeneblast file containing NCBI references to be downloaded.")
def main(mgbfile, outputdir):
    """ Console script for BioCompass"""
    download_hits(mgbfile, outputdir)


if __name__ == "__main__":
    main()
