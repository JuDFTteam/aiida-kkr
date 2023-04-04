#!/usr/bin/env python

from builtins import object
import pytest
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, aiida_profile
from ..conftest import import_with_migration


# tests
def test_parse_kkr_calc(aiida_profile):
    """
    ...
    """
    from aiida.orm import load_node
    from aiida_kkr.parsers.kkr import KkrParser
    import_with_migration('files/db_dump_kkrcalc.tar.gz')
    kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')
    parser = KkrParser(kkr_calc)
    out = parser.parse(debug=False)
    assert out is None
    out_dict = parser.outputs.output_parameters.get_dict()
    # remove this one error message because it is expected
    err = [i for i in out_dict['parser_errors'] if 'OUTPUT_2' not in i]
    assert err == []
