#!/usr/bin/env python

from builtins import object
import pytest
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, aiida_profile
from ..conftest import import_with_migration


# tests
def test_parse_KKRnano_calc(aiida_profile, data_regression):
    """
    Import a KKRnano calculation and parse the outcome
    """
    from aiida.orm import load_node
    from aiida_kkr.parsers.kkrnano import KKRnanoParser

    imported_nodes = import_with_migration('files/db_dump_KKRnano.aiida')
    kkrnano_calc = load_node('da234580-b279-490b-972a-cdf9254c9278')
    parser = KKRnanoParser(kkrnano_calc)
    out = parser.parse(debug=False)
    assert out is None
    out_dict = parser.outputs.output_parameters.get_dict()
    data_regression.check(out_dict['rms_all_iterations'])
