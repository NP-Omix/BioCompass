#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
"""

from click.testing import CliRunner
from BioCompass import cli
def test_importcli():
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    #assert 'BioCompass.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    #assert 'Usage: main [OPTIONS] MGBFILE' in help_result.output


import pandas as pd
from BioCompass.BioCompass import find_category_from_product
def test_find_category_from_product():
    df = pd.DataFrame(
            {'product': ['permease', 'bla bla permease bla bla']})
    out = find_category_from_product(df)
    assert (out['category'] == 'transporter').all()

def test_find_category_from_product_hypothetical():
    """ Uncataloged product should return hypothetical
    """
    df = pd.DataFrame(
            {'product': ['ipsisLitteris']})
    out = find_category_from_product(df)
    assert (out['category'] == 'hypothetical').all()
