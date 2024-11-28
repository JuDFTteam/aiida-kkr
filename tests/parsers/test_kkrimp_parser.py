#!/usr/bin/env python

from builtins import object
import pytest
from aiida import orm
from ..conftest import import_with_migration

# tests


def test_parse_kkrimp_calc(aiida_profile):
    """
    simple Cu noSOC, FP, lmax2
    """
    from aiida_kkr.parsers.kkrimp import KkrimpParser
    group_pk = import_with_migration('data_dir/kkrimp_full_wc.aiida')
    kkrimp_calc = [
        i for i in orm.load_group(group_pk).nodes if i.label == 'KKRimp calculation step 2 (IMIX=0, Zimp: [30.0])'
    ][0]
    print(kkrimp_calc, kkrimp_calc.label)

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
    from aiida_kkr.parsers.kkrimp import KkrimpParser
    import_with_migration('files/export_kkrimp_calc.aiida')
    kkrimp_calc = orm.load_node('b9a2b29a-e250-4992-ae6a-579b733ad1f8')
    parser = KkrimpParser(kkrimp_calc)
    out = parser.parse(debug=False, doscalc=False)
    assert out is None
    out_dict = parser.outputs.output_parameters.get_dict()
    assert out_dict['parser_errors'] == []
