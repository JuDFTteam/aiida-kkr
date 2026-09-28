#!/usr/bin/env python
"""kkr_imp_sub_wc sets the non-spherical cutoffs on the simple-mixing step only when asked (no database needed)."""

from types import SimpleNamespace
from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc


def _stub(factor):
    return SimpleNamespace(ctx=SimpleNamespace(pot_ns_cutoff_factor=factor, convergence_criterion=1e-7))


def test_pot_ns_cutoffs():
    """None leaves KKRimp's defaults; a factor sets both keys on simple mixing and restores the defaults afterwards."""
    cutoffs = kkr_imp_sub_wc._pot_ns_cutoffs
    assert cutoffs(_stub(None), simple_mixing=True) == {}
    assert cutoffs(_stub(None), simple_mixing=False) == {}
    assert cutoffs(_stub(0.1), simple_mixing=True) == {'POT_NS_CUTOFF': 0.1 * 1e-7, 'POT_NS_WRITE_CUTOFF': 0.1 * 1e-7}
    # live test 888237: without these, the Anderson step inherited 1e-8 from its parent calculation
    assert cutoffs(_stub(0.1), simple_mixing=False) == {'POT_NS_CUTOFF': 0.1 * 1e-7, 'POT_NS_WRITE_CUTOFF': 1e-7}
