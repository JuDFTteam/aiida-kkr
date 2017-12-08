# -*- coding: utf-8 -*-
"""
@author: ruess
"""

import pytest
from aiida_kkr.tools.kkrparser_functions import parse_kkr_outputfile


class Test_voronoi_parser_functions():
    """
    Tests for the kkr parser functions
    """
    #some global definitions
    global dref, grouping_ref, outfile, outfile_0init, outfile_000, timing_file, potfile_out
    dref = {}
    grouping_ref = ['energy_contour_group', 'warnings_group', 'ewald_sum_group', 
                    'timings_group', 'core_states_group', 'convergence_group',
                    'kmesh_group', 'symmetries_group', 'code_info_group']
    path0 = '../tests/files/kkr/'
    outfile = path0+'out_kkr'
    outfile_0init = path0+'output.0.txt'
    outfile_000 = path0+'output.000.txt'
    timing_file = path0+'out_timing.000.txt'
    potfile_out = path0+'out_potential'
    
    
    def test_complete_kkr_output(self):
        """
        Parse complete output of kkr calculation
        """
        out_dict = {}
        success, msg_list, out_dict = parse_kkr_outputfile(out_dict, outfile, outfile_0init, outfile_000, timing_file, potfile_out)
        out_dict['parser_warnings'] = msg_list
        assert success
        assert out_dict == dref
        assert msg_list == []
        groups = [i for i in out_dict.keys() if 'group' in i]
        assert set(groups) == set(grouping_ref)

    def test_missing_outfile(self):
        """
        Parse kkr output where out_kkr is missing. Compares error messages and rest of out_dict
        """
        #TODO implement test
        assert 1==2

    def test_missing_outfile0init(self):
        """
        Parse kkr output where output.0.txt is missing. Compares error messages and rest of out_dict
        """
        #TODO implement test
        assert 1==2

    def test_missing_outfile000(self):
        """
        Parse kkr output where output.000.txt is missing. Compares error messages and rest of out_dict
        """
        #TODO implement test
        assert 1==2

    def test_missing_timingfile(self):
        """
        Parse kkr output where out_timing.000.txt is missing. Compares error messages and rest of out_dict
        """
        #TODO implement test
        assert 1==2

    def test_missing_potfile(self):
        """
        Parse kkr output where out_potential is missing. Compares error messages and rest of out_dict
        """
        #TODO implement test
        assert 1==2