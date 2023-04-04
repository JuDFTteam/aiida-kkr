#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import pytest

from builtins import object
from aiida_kkr.tools.tools_kkrimp import modify_potential
from masci_tools.io.common_functions import open_general
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test


class Test_modify_potential(object):
    """ Tests for the modify_potential class functions. """

    def test_neworder_potential_filehandle(self):
        pot1 = 'files/kkr/kkr_run_slab_nosoc/out_potential'
        pot_out1 = 'test_pot1'
        pot_out2 = 'test_pot2'
        pot_out3 = 'test_pot3'
        pot_out4 = 'test_pot4'
        pot1_fh = open_general(pot1)
        pot_out_fh3 = open_general(pot_out3, 'w')
        pot_out_fh4 = open_general(pot_out4, 'w')
        # test using file paths
        modify_potential().neworder_potential(pot1, pot_out1, [0, 0, 0])
        with open(pot_out1) as f:
            po1 = f.readlines()
        # test file handle on input
        modify_potential().neworder_potential(pot1_fh, pot_out2, [0, 0, 0])
        with open(pot_out2) as f:
            po2 = f.readlines()
        # test file handle on output
        modify_potential().neworder_potential(pot1, pot_out_fh3, [0, 0, 0])
        with open(pot_out3) as f:
            po3 = f.readlines()
        # test file handle both
        modify_potential().neworder_potential(pot1_fh, pot_out_fh4, [0, 0, 0])
        with open(pot_out4) as f:
            po4 = f.readlines()
        # verify result
        print(len(po1), len(po2), len(po3), len(po4))
        assert po1 == po2
        assert po1 == po3
        assert po1 == po4

    def test_shapefun_from_scoef(self, file_regression):
        shapefun_path = '../tests/files/mod_pot/test2/shapefun'
        scoefpath = '../tests/files/mod_pot/test2/scoef'
        atom2shapes = [1]
        shapefun_new = '../tests/files/mod_pot/test2/shapefun_new'
        modify_potential().shapefun_from_scoef(scoefpath, shapefun_path, atom2shapes, shapefun_new)
        with open(shapefun_new) as _f:
            txt = ''.join([line for line in _f])
        file_regression.check(txt)

    def test_neworder_potential_no_replace(self, file_regression):
        path = '../tests/files/mod_pot/test1/'
        pot = path + 'pot'
        out_pot = path + 'pot_new'
        neworder = [0, 1, 2]
        # test 1: neworder_potential standard
        modify_potential().neworder_potential(pot, out_pot, neworder)
        with open(out_pot) as _f:
            txt = ''.join([line for line in _f])
        file_regression.check(txt)

    def test_neworder_potential_with_replace(self, file_regression):
        path = '../tests/files/mod_pot/test1/'
        pot1 = path + 'out_potential'
        pot2 = path + 'pot'
        out_pot = path + 'pot_new'
        neworder = [0, 1, 2]
        # test 2: neworder_potential with replace from second potential
        replace_newpos = [[0, 0], [2, 0]]
        modify_potential().neworder_potential(pot1, out_pot, neworder, potfile_2=pot2, replace_from_pot2=replace_newpos)
        # now compare to reference:
        with open(out_pot) as _f:
            txt = ''.join([line for line in _f])
        file_regression.check(txt)
