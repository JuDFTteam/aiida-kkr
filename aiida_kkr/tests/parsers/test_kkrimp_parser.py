#!/usr/bin/env python

from __future__ import absolute_import
from builtins import object
import pytest

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkrimp_parser(object):
    """
    Tests for the kkrimp calculation
    """

    def test_parse_kkrimp_calc(self):
        """
        simple Cu noSOC, FP, lmax2
        """
        from aiida.orm import load_node
        from aiida_kkr.parsers.kkrimp import KkrimpParser
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrimp_out.tar.gz')
        kkrimp_calc = load_node('eab8db1b-2cc7-4b85-a524-0df4ff2b7da6')
        parser = KkrimpParser(kkrimp_calc)
        parser.parse()


if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   t = Test_kkrimp_parser()
   t.test_parse_kkrimp_calc()
