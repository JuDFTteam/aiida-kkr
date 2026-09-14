#!/usr/bin/env python
"""Routing of a diverged KKRimp calculation in kkr_imp_sub_wc (no database needed)."""

from types import SimpleNamespace
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc


def _workchain_stub(exit_status):
    """Minimal stand-in for a kkr_imp_sub_wc instance whose last KKRimp calculation failed."""
    output_parameters = SimpleNamespace(get_dict=lambda: {'nonfinite_values': ['convergence_group.rms']})
    last_calc = SimpleNamespace(
        is_finished_ok=False, exit_status=exit_status, outputs=SimpleNamespace(output_parameters=output_parameters)
    )
    return SimpleNamespace(
        ctx=SimpleNamespace(calcs=[], last_calc=last_calc),
        exit_codes=kkr_imp_sub_wc.exit_codes,
        report=lambda message: None,
    )


def test_inspect_kkrimp_nonfinite_output():
    """A calculation that stopped with ERROR_NONFINITE_OUTPUT (303) ends the workchain with 134."""
    stub = _workchain_stub(303)
    assert kkr_imp_sub_wc.inspect_kkrimp(stub).status == 134
    assert stub.ctx.kkrimp_step_success is False


def test_inspect_kkrimp_other_failure():
    """Any other failed calculation still ends the workchain with 130."""
    assert kkr_imp_sub_wc.inspect_kkrimp(_workchain_stub(302)).status == 130
