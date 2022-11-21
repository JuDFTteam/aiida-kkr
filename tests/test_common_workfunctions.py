# -*- coding: utf-8 -*-
"""
@author: ruess
"""

from builtins import object
import pytest
from numpy import sort
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test
from .conftest import import_with_migration


@pytest.mark.usefixtures('aiida_profile')
class Test_common_workfunctions(object):
    """
    Tests for the common workfunctions from tools.common_workfunctions,
    i.e. functions commonly used in this plugin that depend on aiida stuff to work
    """

    def test_generate_inputcard_from_structure(self):
        from aiida_kkr.tools import generate_inputcard_from_structure
        from aiida.orm import StructureData, Dict

        s = StructureData(cell=[[0.5, 0.5, 0], [1, 0, 0], [0, 0, 1]])
        s.append_atom(position=[0, 0, 0], symbols='Fe')
        p = Dict({'LMAX': 2, 'NSPIN': 2, 'RMAX': 10, 'GMAX': 100})
        generate_inputcard_from_structure(p, s, 'inputcard')
        fhandle = open('inputcard2', 'w')
        generate_inputcard_from_structure(p, s, fhandle)
        txt = open('inputcard', 'r').readlines()  # from path
        txt2 = open('inputcard2', 'r').readlines()  # from fhandle
        ref = [
            'ALATBASIS= 1.889726124626\n', 'BRAVAIS\n', '0.500000000000 0.500000000000 0.000000000000\n',
            '1.000000000000 0.000000000000 0.000000000000\n', '0.000000000000 0.000000000000 1.000000000000\n',
            'NAEZ= 1\n', '<RBASIS>\n', '0.000000000000 0.000000000000 0.000000000000\n', 'CARTESIAN= True\n',
            '<ZATOM>\n', '26.000000000000\n', 'NSPIN= 2\n', 'LMAX= 2\n', 'RMAX=     10.000000000000\n',
            'GMAX= 100.000000000000\n'
        ]
        done = False
        while not done:
            try:
                txt.remove('\n')
            except ValueError:
                done = True
        assert len(txt) == len(ref)
        txt.sort()
        ref.sort()
        print(txt, ref)
        for i in range(len(txt)):
            print(i, txt[i], ref[i])
            assert set(txt[i].split()) == set(ref[i].split())
        # check if fhandle output also works
        done = False
        while not done:
            try:
                txt2.remove('\n')
            except ValueError:
                done = True
        assert len(txt2) == len(ref)
        txt2.sort()
        for i in range(len(txt)):
            assert set(txt2[i].split()) == set(ref[i].split())

    def test_check_2Dinput_consistency_1(self):
        # case 1: 3D structure and no 2D params input
        from aiida_kkr.tools import check_2Dinput_consistency
        from aiida.orm import StructureData, Dict

        s = StructureData(cell=[[0.5, 0.5, 0], [1, 0, 0], [0, 0, 1]])
        s.append_atom(position=[0, 0, 0], symbols='Fe')
        p = Dict({'INTERFACE': False})
        input_check = check_2Dinput_consistency(s, p)
        assert input_check[0]
        assert input_check[1] == '2D consistency check complete'

    def test_check_2Dinput_consistency_2(self):
        # case 2: 2D structure and 2D params input
        from aiida_kkr.tools.common_workfunctions import check_2Dinput_consistency
        from aiida.orm import StructureData, Dict

        s = StructureData(cell=[[0.5, 0.5, 0], [1, 0, 0], [0, 0, 1]])
        s.append_atom(position=[0, 0, 0], symbols='Fe')
        s.set_pbc((True, True, False))
        p = Dict({
            'INTERFACE': True,
            '<NRBASIS>': 1,
            '<RBLEFT>': [0, 0, 0],
            '<RBRIGHT>': [0, 0, 0],
            'ZPERIODL': [0, 0, 0],
            'ZPERIODR': [0, 0, 0],
            '<NLBASIS>': 1
        })
        input_check = check_2Dinput_consistency(s, p)
        assert input_check[0]
        assert input_check[1] == '2D consistency check complete'

    def test_check_2Dinput_consistency_3(self):
        # case 3: 2D structure but incomplete 2D input parameters given
        from aiida_kkr.tools import check_2Dinput_consistency
        from aiida.orm import StructureData, Dict

        s = StructureData(cell=[[0.5, 0.5, 0], [1, 0, 0], [0, 0, 1]])
        s.append_atom(position=[0, 0, 0], symbols='Fe')
        s.set_pbc((True, True, False))
        p = Dict({
            'INTERFACE': True,
            '<NRBASIS>': 1,
        })
        input_check = check_2Dinput_consistency(s, p)
        assert not input_check[0]
        assert list(input_check[1]).sort() == list(
            '2D info given in parameters but structure is 3D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}'
            .format(True, False, ['INTERFACE', '<NRBASIS>'])
        ).sort()

    def test_check_2Dinput_consistency_4(self):
        # case 3: 2D structure but interface parameter set to False
        from aiida_kkr.tools import check_2Dinput_consistency
        from aiida.orm import StructureData, Dict

        s = StructureData(cell=[[0.5, 0.5, 0], [1, 0, 0], [0, 0, 1]])
        s.append_atom(position=[0, 0, 0], symbols='Fe')
        s.set_pbc((True, True, False))
        p = Dict({
            'INTERFACE': False,
            '<NRBASIS>': 1,
            '<RBLEFT>': [0, 0, 0],
            '<RBRIGHT>': [0, 0, 0],
            'ZPERIODL': [0, 0, 0],
            'ZPERIODR': [0, 0, 0],
            '<NLBASIS>': 1
        })
        input_check = check_2Dinput_consistency(s, p)
        assert not input_check[0]
        assert input_check[1] == "'INTERFACE' parameter set to False but structure is 2D"

    def test_check_2Dinput_consistency_5(self):
        # case 5: 3D structure but 2D params given
        from aiida_kkr.tools import check_2Dinput_consistency
        from aiida.orm import StructureData, Dict

        s = StructureData(cell=[[0.5, 0.5, 0], [1, 0, 0], [0, 0, 1]])
        s.append_atom(position=[0, 0, 0], symbols='Fe')
        s.set_pbc((True, True, True))
        p = Dict({
            'INTERFACE': True,
            '<NRBASIS>': 1,
            '<RBLEFT>': [0, 0, 0],
            '<RBRIGHT>': [0, 0, 0],
            'ZPERIODL': [0, 0, 0],
            'ZPERIODR': [0, 0, 0],
            '<NLBASIS>': 1
        })
        input_check = check_2Dinput_consistency(s, p)
        assert not input_check[0]
        assert list(input_check[1]).sort() == list(
            '3D info given in parameters but structure is 2D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}'
            .format(
                False, True, ['ZPERIODL', '<NRBASIS>', '<RBLEFT>', 'INTERFACE', '<NLBASIS>', 'ZPERIODR', '<RBRIGHT>']
            )
        ).sort()

    def test_vca_check(self):
        from aiida_kkr.tools.common_workfunctions import vca_check
        pass

    def test_kick_out_corestates_wf(self):
        from aiida_kkr.tools import kick_out_corestates_wf
        from aiida.orm import load_node, Float
        from masci_tools.io.parsers.kkrparser_functions import get_corestates_from_potential
        import numpy as np

        import_with_migration('files/kick_out_corestates_input.tar.gz')
        sfd_potential_node = load_node('933ddebb-e72f-43b0-aca5-cd0a109da75f')
        emin_node = Float(-1.2)
        pot_removed_core = kick_out_corestates_wf(sfd_potential_node, emin_node)

        # compare list of core states before and after kick_out calcfunction
        with sfd_potential_node.open(sfd_potential_node.filename) as potfile:
            ncore0, ecore0, lcore0 = get_corestates_from_potential(potfile)
        with pot_removed_core.open(pot_removed_core.filename) as potfile:
            ncore1, ecore1, lcore1 = get_corestates_from_potential(potfile)

        # compare number of core state
        print(np.array(ncore0) - np.array(ncore1))
        assert np.sum(np.array(ncore0) - np.array(ncore1)) == 2  # check if exactly two cores states have been removed

    """
    def test_prepare_VCA_structure_wf(self):
        #TODO: implement check
        from aiida_kkr.tools.common_workfunctions import prepare_VCA_structure_wf
        pass


    def test_prepare_2Dcalc_wf(self):
        #TODO: implement check
        from aiida_kkr.tools.common_workfunctions import prepare_2Dcalc_wf
        pass


    def test_test_and_get_codenode(self):
        #TODO: implement check
        from aiida_kkr.tools.common_workfunctions import test_and_get_codenode
        pass


    def test_get_inputs_kkr(self):
        #TODO: implement check
        from aiida_kkr.tools.common_workfunctions import get_inputs_kkr
        assert 1==2


    def test_get_inputs_voronoi(self):
        #TODO: implement check
        from aiida_kkr.tools.common_workfunctions import get_inputs_voronoi
        assert 1==2


    def test_get_parent_paranode(self):
        #TODO: implement check
        from aiida_kkr.tools.common_workfunctions import get_parent_paranode
        assert 1==2
    """


#"""
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    from aiida.orm import StructureData, Dict

    t = Test_common_workfunctions()

    t1 = t.test_update_params_wf()
    #t2 = t.test_prepare_VCA_structure_wf()
    #t3 = t.test_prepare_2Dcalc_wf()
    #t4 = t.test_test_and_get_codenode()
    #t5 = t.test_get_inputs_kkr()
    #t7 = t.test_get_parent_paranode()
#"""
