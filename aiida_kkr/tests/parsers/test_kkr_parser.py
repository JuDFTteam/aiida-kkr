#!/usr/bin/env python

from __future__ import absolute_import
from builtins import object
import pytest

from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test

# tests
class Test_kkr_parser(object):
    """
    Tests for the kkr parser
    """

    def test_parse_kkr_calc(self, clear_database):
        """
        ...
        """
        from aiida.orm import load_node
        from aiida_kkr.parsers.kkr import KkrParser
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')
        parser = KkrParser(kkr_calc)
        parser.parse()
