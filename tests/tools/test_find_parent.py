#!/usr/bin/env python

import pytest
from aiida.orm import load_node
from aiida_kkr.tools.find_parent import find_parent_structure, get_calc_from_remote
from ..conftest import import_with_migration


def test_find_parent_structure():
    """
    find parent structure and Voronoi parent calculation from KKR calculation
    """

    # load necessary files from db_dump files
    import_with_migration('files/db_dump_kkrcalc.tar.gz')
    kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

    # now find voronoi parent and structure
    struc, voro_parent = find_parent_structure(kkr_calc)

    # check result
    assert struc.uuid == 'e51ee6a1-bd27-4901-9612-7bac256bf117'
    assert voro_parent.uuid == '559b9d9b-3525-402e-9b24-ecd8b801853c'


def test_find_structure_kkrimp():
    """
    find parent structure from a KkrimpCalculation
    """
    import_with_migration('files/db_dump_kkrimp_out.tar.gz')
    kkrimp_calc = load_node('eab8db1b-2cc7-4b85-a524-0df4ff2b7da6')

    # now find voronoi parent and structure
    struc, voro_parent = find_parent_structure(kkrimp_calc)

    # check result
    assert struc.uuid == 'e51ee6a1-bd27-4901-9612-7bac256bf117'
    assert voro_parent.uuid == '559b9d9b-3525-402e-9b24-ecd8b801853c'


def test_get_calc_from_remote():
    """
    find parent calc from remote
    """

    # load necessary files from db_dump files
    import_with_migration('files/db_dump_kkrcalc.tar.gz')
    kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

    # now find voronoi parent and structure
    parent_from_remote = get_calc_from_remote(kkr_calc.outputs.remote_folder)

    # check result
    assert parent_from_remote.uuid == kkr_calc.uuid
