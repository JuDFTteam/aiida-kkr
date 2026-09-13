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
    # version information must survive parsing, it is the only way to tell which parser produced a node
    assert out_dict['parser_version'] == parser._ParserVersion
    assert 'calculation_plugin_version' in out_dict
    assert 'nonfinite_values' not in out_dict


def test_parse_kkr_nonfinite(aiida_profile):
    """
    diverged KKR host calculation whose output contains NaN values

    Retrieved files of a real run (aiida-core 1.5.2, aiida-kkr 1.1.11-dev4, step 5 of a kkr_scf_wc with
    Broyden mixing) that blows up in its 3rd and last iteration. masci-tools parses it with success, so the
    non-finite values have to be caught independently of the parser's verdict.
    """
    from aiida.orm import load_node
    from aiida_kkr.parsers.kkr import KkrParser
    import_with_migration('files/kkr/nonfinite_kkrhost_2128.aiida')
    kkr_calc = load_node('29d772ea-88dc-4c19-aac9-ceb65f569ca6')
    parser = KkrParser(kkr_calc)
    exit_code = parser.parse(debug=False)
    assert exit_code.status == 304

    # storing is where the unsanitized output raised 'nan and inf/-inf can not be serialized to the database'
    out_dict = parser.outputs.output_parameters.store().get_dict()

    assert out_dict['parser_version'] == parser._ParserVersion
    assert len(out_dict['nonfinite_values']) == 15
    for path in (
        'energy', 'total_energy_Ry', 'convergence_group.rms',
        'convergence_group.spin_moment_per_atom_all_iterations[3][0]'
    ):
        assert path in out_dict['nonfinite_values']
    convergence = out_dict['convergence_group']
    assert convergence['first_nonfinite_iteration_index'] == 2
    assert len(convergence['rms_all_iterations']) == 3
    assert convergence['rms_all_iterations'][2] is None
    assert out_dict['energy'] is None and out_dict['total_energy_Ry'] is None
    assert any('nonfinite_values' in warning for warning in out_dict['parser_warnings'])
