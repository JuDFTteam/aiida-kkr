#!/usr/bin/env python
"""Tests for the non-finite sanitizer used by the KKR and KKRimp parsers (no database needed)."""

import math
import numpy as np
from aiida.orm.implementation.utils import clean_value
from aiida_kkr.tools.dict_util import sanitize_nonfinite, record_nonfinite

NAN, INF = float('nan'), float('inf')


def _ragged_spin_moments():
    """Shape of convergence_group.total_spin_moment_all_iterations in a real diverged KKRimp run
    (Ba:Cu, calc 871304): entries [0] and [1] are nested three deep, [2] is a bare scalar."""
    block1 = [[0.0, 0.0, 0.0] for _ in range(513)] + [[0.0, 0.0, INF]]
    return [[[0.1, 0.2, NAN]], block1, NAN]


def _out_dict():
    return {
        'energy': NAN,
        'total_charge_per_atom': [1.0, -INF],
        'number_of_atoms_in_unit_cell': 1,
        'calculation_converged': False,
        'convergence_group': {
            'rms': NAN,
            'rms_all_iterations': [0.3, 0.2, 0.31, NAN],
            'total_spin_moment_all_iterations': _ragged_spin_moments(),
            'nested': {
                'tuple': (np.float64(NAN), 2.0)
            },
        },
        'array': np.array([1.0, INF]),
        'parser_errors': [],
    }


def test_sanitize_nonfinite_paths_and_replacement():
    """Every non-finite float is found at any depth, replaced by None in place, lists keep their length."""
    out = _out_dict()
    paths = sanitize_nonfinite(out)

    assert paths == [
        'energy',
        'total_charge_per_atom[1]',
        'convergence_group.rms',
        'convergence_group.rms_all_iterations[3]',
        'convergence_group.total_spin_moment_all_iterations[0][0][2]',
        'convergence_group.total_spin_moment_all_iterations[1][513][2]',
        'convergence_group.total_spin_moment_all_iterations[2]',
        'convergence_group.nested.tuple[0]',
        'array[1]',
    ]
    conv = out['convergence_group']
    assert conv['rms_all_iterations'] == [0.3, 0.2, 0.31, None]
    assert conv['total_spin_moment_all_iterations'][0] == [[0.1, 0.2, None]]
    assert len(conv['total_spin_moment_all_iterations'][1]) == 514
    assert conv['total_spin_moment_all_iterations'][1][513] == [0.0, 0.0, None]
    assert conv['total_spin_moment_all_iterations'][2] is None
    assert conv['nested']['tuple'] == (None, 2.0)
    assert out['array'] == [1.0, None]
    # integers and bools are untouched
    assert out['number_of_atoms_in_unit_cell'] == 1 and out['calculation_converged'] is False
    # what AiiDA does when storing the Dict must now pass
    clean_value(out)


def test_sanitize_nonfinite_clean_input():
    """A dict without non-finite values is returned unchanged and nothing is recorded."""
    out = {'rms': 1e-8, 'rms_all_iterations': [1.0, 1e-8], 'parser_errors': []}
    paths = sanitize_nonfinite(out)
    record_nonfinite(out, paths)
    assert paths == []
    assert out == {'rms': 1e-8, 'rms_all_iterations': [1.0, 1e-8], 'parser_errors': []}


def test_record_nonfinite():
    """The replaced paths, the first diverged iteration and a warning end up in the output."""
    out = _out_dict()
    out['convergence_group']['rms_all_iterations'] = [0.3, NAN, 0.2, NAN]
    paths = sanitize_nonfinite(out)
    record_nonfinite(out, paths)

    assert out['nonfinite_values'] == paths
    assert out['convergence_group']['first_nonfinite_iteration_index'] == 1
    assert len(out['parser_warnings']) == 1 and 'nonfinite_values' in out['parser_warnings'][0]
    assert out['parser_errors'] == []
    clean_value(out)


def test_record_nonfinite_without_rms_list():
    """No iteration index is invented when rms_all_iterations is finite or absent."""
    out = {'energy': NAN, 'convergence_group': {'rms_all_iterations': [0.3, 0.2]}, 'parser_warnings': ['old']}
    record_nonfinite(out, sanitize_nonfinite(out))
    assert out['nonfinite_values'] == ['energy']
    assert 'first_nonfinite_iteration_index' not in out['convergence_group']
    assert out['parser_warnings'][0] == 'old' and len(out['parser_warnings']) == 2
    assert not any(isinstance(v, float) and math.isnan(v) for v in out.values())
