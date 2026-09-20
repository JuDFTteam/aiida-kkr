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


def get_test_struc_open():
    """
    Open bcc lattice, a = 5.028 Ang, single site at the origin.

    Its 79-atom cluster does not fit into the 10 Ang default search radius.
    """
    alat = 5.028
    s = StructureData(
        cell=[[-alat / 2, alat / 2, alat / 2], [alat / 2, -alat / 2, alat / 2], [alat / 2, alat / 2, -alat / 2]]
    )
    s.append_atom(position=[0, 0, 0], symbols='Fe')
    return s


def get_test_struc_two_site(shift=(0., 0., 0.)):
    """
    Two-site cubic cell, a = 4 Ang, Cu at the origin and Au at the body center.

    `shift` translates both sites by the same vector, which moves them away from the
    origin without changing the physics.
    """
    shift = np.array(shift, dtype=float)
    s = StructureData(cell=[[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]])
    s.append_atom(position=list(shift), symbols='Cu')
    s.append_atom(position=list(shift + 2.), symbols='Au')
    return s


def test_find_cluster_radius_grows_search_radius():
    """
    The search radius has to grow when the requested cluster does not fit into it.

    This lattice has only 58 neighbors within the 10 Ang default, so asking for 79 atoms
    used to index past the end of the neighbor list (issue #182).
    """
    r0, ncls = find_cluster_radius(get_test_struc_open(), 79)
    assert 10. < r0 < 12.
    assert min(ncls) >= 79


def test_find_cluster_radius_off_origin():
    """
    Neighbor distances are measured from the central site, not from the origin.

    Two consequences of measuring from the origin, on the same two-site cell: the answer
    changed when all sites were translated, and the returned radius held fewer atoms than
    were asked for (issue #182).
    """
    r0_origin = find_cluster_radius(get_test_struc_two_site(), 15)[0]
    r0_shifted = find_cluster_radius(get_test_struc_two_site(shift=(1., 0., 0.)), 15)[0]
    assert np.round(r0_origin, 5) == np.round(r0_shifted, 5) == 4.

    ncls = find_cluster_radius(get_test_struc_two_site(shift=(1., 0., 0.)), 18)[1]
    assert min(ncls) >= 18


def test_find_cluster_radius_nclsmin_too_small():
    """
    A cluster of less than two atoms has no neighbor distance to return.
    """
    with pytest.raises(ValueError):
        find_cluster_radius(get_test_struc(), 1)


def test_find_cluster_radius_exactly_enough_neighbors():
    """
    A starting radius holding exactly the needed neighbors is used as it is.

    One neighbor short of that, the radius has to grow rather than index past the end of
    the list. The old code got the first case right only by accident, by reading the last
    element of the array.
    """
    s = get_test_struc()
    # this radius holds exactly 18 neighbors of this cell, i.e. a cluster of 19 atoms
    r18 = 1.2247449
    assert len(s.get_pymatgen().get_all_neighbors(r18)[0]) == 18

    r0, ncls = find_cluster_radius(s, 19, Rclsmax=r18)
    assert np.round(r0, 6) == 1.224745
    assert ncls[0] == 19

    r0, ncls = find_cluster_radius(s, 20, Rclsmax=r18)
    assert r0 > r18
    assert ncls[0] >= 20
