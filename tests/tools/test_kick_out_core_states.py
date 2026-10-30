#!/usr/bin/env python
# coding: utf-8
"""
Test for kick_out_corestates: core states with l > LMAX are kept
"""

from pathlib import Path
import pytest
from masci_tools.io.common_functions import get_corestates_from_potential
from aiida_kkr.tools.kick_out_core_states import kick_out_corestates

# first block (atom 1, spin down) of the input potential of a KKR host scf step for compressed Hf, 12 core states,
# highest are 5p (l=1, -1.52 Ry) and 4f (l=3, -0.07 Ry)
POTFILE = Path(__file__).parent.parent / 'files' / 'kick_out_corestates' / 'potential_Hf_0.8_kkr98841'
EMIN = -0.3  # EMIN of that run, only the 4f lies inside the contour


@pytest.mark.parametrize('lmax, n_removed', [(2, 0), (3, 1), (None, 1)])
def test_kick_out_corestates_lmax(tmp_path, lmax, n_removed):
    potfile_out = tmp_path / 'potential'
    num_deleted = kick_out_corestates(str(POTFILE), str(potfile_out), EMIN, lmax=lmax)

    n0, _, _ = get_corestates_from_potential(str(POTFILE))
    assert num_deleted == n_removed * len(n0)
    if n_removed == 0:
        # nothing written if no core state is removed
        assert not potfile_out.exists()
    else:
        n1, e1, l1 = get_corestates_from_potential(str(potfile_out))
        assert list(n1) == [n - n_removed for n in n0]
        for ll, ee in zip(l1, e1):
            assert 3 not in ll
            assert (ee <= EMIN).all()
