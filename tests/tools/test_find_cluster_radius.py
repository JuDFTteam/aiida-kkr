#!/usr/bin/env python
# coding: utf-8
"""
Test for find_cluster_radius functionality
"""

import pytest
import numpy as np
from aiida.orm import StructureData
from aiida_kkr.tools.find_cluster_radius import find_cluster_radius_old, find_cluster_radius


def get_test_struc():
    """
    Test structure
    """
    s = StructureData(cell=[[0.5, 0.5, 0], [1, 0, 0], [0, 0, 1]])
    s.append_atom(position=[0, 0, 0], symbols='Fe')
    return s


def test_find_cluster_radius_old():
    """
    Test for old version
    """
    s = get_test_struc()
    r, r_alat = find_cluster_radius_old(s, 15, n_max_box=50, nbins=100)

    saux = StructureData(cell=s.cell)
    for site in s.sites:
        kind = s.get_kind(site.kind_name)
        saux.append_atom(position=site.position, symbols=kind.symbols, weights=kind.weights)
    saux.pbc = (True, True, True)

    ps = saux.get_pymatgen()

    neighbors_all = ps.get_all_neighbors(r)
    assert r > 1.285
    assert len(neighbors_all[0]) == 18


def test_find_cluster_radius():
    """
    Test for new version
    """
    r0, ncls = find_cluster_radius(get_test_struc(), 15)
    assert np.round(r0, 3) == 1.225
    assert ncls[0] == 19
