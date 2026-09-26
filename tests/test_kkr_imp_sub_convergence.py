#!/usr/bin/env python
"""kkr_imp_sub_wc counts a run as converged only at convergence_criterion (no database needed)."""

from types import SimpleNamespace
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc


def _workchain_stub(qbound, **ctx):
    """Minimal stand-in for a kkr_imp_sub_wc instance whose last KKRimp calculation ran with QBOUND=qbound."""
    parameters = SimpleNamespace(get_dict=lambda: {'QBOUND': qbound})
    last_calc = SimpleNamespace(is_finished_ok=True, inputs=SimpleNamespace(parameters=parameters))
    stub = SimpleNamespace(
        ctx=SimpleNamespace(
            last_calc=last_calc, convergence_criterion=1e-7, loop_count=2, max_number_runs=5, successful=False, **ctx
        ),
        exit_codes=kkr_imp_sub_wc.exit_codes,
        report=lambda message: None,
    )
    stub._last_calc_qbound = lambda: kkr_imp_sub_wc._last_calc_qbound(stub)
    return stub


def _output(calculation_converged):
    return {'convergence_group': {'calculation_converged': calculation_converged}}


def test_reached_convergence_criterion():
    """Reaching the simple-mixing QBOUND is not convergence; reaching convergence_criterion is."""
    reached = kkr_imp_sub_wc._reached_convergence_criterion
    assert not reached(_workchain_stub(1e-2), _output(True))
    assert reached(_workchain_stub(1e-7), _output(True))
    assert not reached(_workchain_stub(1e-7), _output(False))
    assert not reached(_workchain_stub(None), _output(True))
    assert reached(_workchain_stub(1e-2), {'convergence_group': {'doscalc': True}})


def test_condition_simple_mixing_restart_keeps_iterating():
    """Workchain 878357: after higher accuracy, a 1-step simple-mixing calc at QBOUND 1e3 must not end the run."""
    stub = _workchain_stub(1e3, kkr_converged=True, kkr_higher_accuracy=True, kkr_converged_to_criterion=False)
    assert kkr_imp_sub_wc.condition(stub) is True
    assert stub.ctx.successful is False


def test_condition_converged_to_criterion_stops():
    """Reaching convergence_criterion after higher accuracy still ends the run with success."""
    stub = _workchain_stub(1e-7, kkr_converged=True, kkr_higher_accuracy=True, kkr_converged_to_criterion=True)
    assert kkr_imp_sub_wc.condition(stub) is False
    assert stub.ctx.successful is True
