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
