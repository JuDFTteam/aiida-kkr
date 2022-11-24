#!/usr/bin/env python

from builtins import object
import pytest
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, aiida_profile
from ..conftest import import_with_migration

# tests


def test_parse_kkrimp_calc(aiida_profile):
    """
    simple Cu noSOC, FP, lmax2
    """
    from aiida.orm import load_node
    from aiida_kkr.parsers.kkrimp import KkrimpParser
    import_with_migration('files/db_dump_kkrimp_out.tar.gz')
    kkrimp_calc = load_node('eab8db1b-2cc7-4b85-a524-0df4ff2b7da6')
    print(kkrimp_calc.outputs.retrieved.list_object_names())
    parser = KkrimpParser(kkrimp_calc)
    out = parser.parse(debug=False)
    print(out)
    assert out is None
    out_dict = parser.outputs.output_parameters.get_dict()
    assert out_dict['parser_errors'] == []


def test_parse_kkrimp_calc_complex(aiida_profile):
    """
    complex magnetic impurity with SOC
    """
    from aiida.orm import load_node
    from aiida_kkr.parsers.kkrimp import KkrimpParser
    import_with_migration('files/export_kkrimp_calc.tar.gz')
    kkrimp_calc = load_node('7547303b-69b7-4380-b0c0-7440e6c4f2a1')
    parser = KkrimpParser(kkrimp_calc)
    out = parser.parse(debug=False)
    assert out is None
    out_dict = parser.outputs.output_parameters.get_dict()
    assert out_dict['parser_errors'] == []
