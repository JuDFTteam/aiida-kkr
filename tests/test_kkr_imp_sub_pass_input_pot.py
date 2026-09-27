#!/usr/bin/env python
"""kkr_imp_sub_wc hands on the input potential after a one-iteration simple-mixing step (no database needed)."""

from types import SimpleNamespace
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc


def _stub(**params):
    parameters = SimpleNamespace(get_dict=lambda: params)
    return SimpleNamespace(
        ctx=SimpleNamespace(last_calc=SimpleNamespace(inputs=SimpleNamespace(parameters=parameters)))
    )


def _output(converged, iterations):
    return {'convergence_group': {'calculation_converged': converged, 'number_of_iterations': iterations}}


def test_quick_simple_mixing():
    """Only a converged one-iteration IMIX 0 step without a magnetic field qualifies (Mn:Cu oracle, calc 881118)."""
    quick = kkr_imp_sub_wc._quick_simple_mixing
    assert quick(_stub(IMIX=0, QBOUND=1e-2), _output(True, 1))
    assert quick(_stub(QBOUND=1e-2), _output(True, 1))  # IMIX unset means simple mixing
    assert quick(_stub(IMIX=0, HFIELD=[0.0, 0]), _output(True, 1))
    assert not quick(_stub(IMIX=0, HFIELD=[0.02, 5]), _output(True, 1))  # mag_init step: keep its moment
    assert not quick(_stub(IMIX=0), _output(True, 6))
    assert not quick(_stub(IMIX=0), _output(False, 1))
    assert not quick(_stub(IMIX=5), _output(True, 1))
