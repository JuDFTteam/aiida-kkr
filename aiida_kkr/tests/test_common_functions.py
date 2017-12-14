# -*- coding: utf-8 -*-
"""
@author: ruess
"""

import pytest
from aiida_kkr.tools.common_functions import interpolate_dos


class Test_common_functions():
    """
    Tests for the common functions from tools.common_functions
    """

    #test interpol:
    def test_interpolate_dos(self):
        from numpy import load
        d0 = '../tests/files/interpol/' 
        ef, dos, dos_int = interpolate_dos(d0, return_original=True)
        assert ef == 0.5256
        assert (dos == load(d0+'/ref_dos.npy')).all()
        assert (dos_int == load(d0+'/ref_dos_int.npy')).all()
        
    def test_get_alat_from_bravais(self):
        #TODO needs implementing
        assert 1==2

    def test_search_string(self):
        #TODO needs implementing
        assert 1==2

    def test_angles_to_vec(self):
        #TODO needs implementing
        assert 1==2

    def test_vec_to_angles(self):
        #TODO needs implementing
        assert 1==2

    def test_get_version_info(self):
        #TODO needs implementing
        assert 1==2

    def test_get_corestates_from_potential(self):
        #TODO needs implementing
        assert 1==2

    def test_get_highest_core_state(self):
        #TODO needs implementing
        assert 1==2

    def test_generate_inputcard_from_structure(self):
        #TODO needs implementing
        assert 1==2

    def test_check_2Dinput(self):
        #TODO needs implementing
        assert 1==2