# -*- coding: utf-8 -*-
"""
@author: ruess
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import object
import pytest
from six.moves import range
from numpy import sort

@pytest.mark.usefixtures("aiida_env")
class Test_common_workfunctions(object):
    """
    Tests for the common workfunctions from tools.common_workfunctions,
    i.e. functions commonly used in this plugin that depend on aiida stuff to work
    """

    def test_generate_inputcard_from_structure(self):
        from aiida_kkr.tools.common_workfunctions import generate_inputcard_from_structure
        from aiida.plugins import DataFactory
        StructureData = DataFactory('structure')
        Dict = DataFactory('dict')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        p = Dict(dict={'LMAX':2, 'NSPIN':2, 'RMAX':10, 'GMAX':100})
        generate_inputcard_from_structure(p, s, 'inputcard')
        fhandle = open('inputcard2', 'w')
        generate_inputcard_from_structure(p, s, fhandle)
        txt = open('inputcard', 'r').readlines() # from path
        txt2 = open('inputcard2', 'r').readlines() # from fhandle
        ref = ['ALATBASIS= 1.88972612545783\n',
               'BRAVAIS\n',
               '0.50000000000000 0.50000000000000 0.00000000000000\n',
               '1.00000000000000 0.00000000000000 0.00000000000000\n',
               '0.00000000000000 0.00000000000000 1.00000000000000\n',
               'NAEZ= 1\n',
               '<RBASIS>\n',
               '0.00000000000000 0.00000000000000 0.00000000000000\n',
               'CARTESIAN= True\n',
               '<ZATOM>\n',
               '26.00000000000000\n',
               'NSPIN= 2\n',
               'LMAX= 2\n',
               'RMAX=     10.00000000000000\n',
               'GMAX= 100.00000000000000\n']
        done=False
        while not done:
            try:
                txt.remove('\n')
            except ValueError:
                done = True
        assert len(txt)==len(ref)
        txt.sort()
        ref.sort()
        print(txt, ref)
        for i in range(len(txt)):
            print(i, txt[i], ref[i])
            assert set(txt[i].split())==set(ref[i].split())
        # check if fhandle output also works
        done=False
        while not done:
            try:
                txt2.remove('\n')
            except ValueError:
                done = True
        assert len(txt2)==len(ref)
        txt2.sort()
        for i in range(len(txt)):
            assert set(txt2[i].split())==set(ref[i].split())


    def test_check_2Dinput_consistency_1(self):
        # case 1: 3D structure and no 2D params input
        from aiida_kkr.tools.common_workfunctions import check_2Dinput_consistency
        from aiida.plugins import DataFactory
        StructureData = DataFactory('structure')
        Dict = DataFactory('dict')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        p = Dict(dict={'INTERFACE':False})
        input_check = check_2Dinput_consistency(s, p)
        assert input_check[0]
        assert input_check[1] == '2D consistency check complete'

    def test_check_2Dinput_consistency_2(self):
        # case 2: 2D structure and 2D params input
        from aiida_kkr.tools.common_workfunctions import check_2Dinput_consistency
        from aiida.plugins import DataFactory
        StructureData = DataFactory('structure')
        Dict = DataFactory('dict')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        s.set_pbc((True, True, False))
        p = Dict(dict={'INTERFACE':True, '<NRBASIS>':1, '<RBLEFT>':[0,0,0], '<RBRIGHT>':[0,0,0], 'ZPERIODL':[0,0,0], 'ZPERIODR':[0,0,0], '<NLBASIS>':1})
        input_check = check_2Dinput_consistency(s, p)
        assert input_check[0]
        assert input_check[1] == "2D consistency check complete"

    def test_check_2Dinput_consistency_3(self):
        # case 3: 2D structure but incomplete 2D input parameters given
        from aiida_kkr.tools.common_workfunctions import check_2Dinput_consistency
        from aiida.plugins import DataFactory
        StructureData = DataFactory('structure')
        Dict = DataFactory('dict')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        s.set_pbc((True, True, False))
        p = Dict(dict={'INTERFACE':True, '<NRBASIS>':1,})
        input_check = check_2Dinput_consistency(s, p)
        assert not input_check[0]
        assert input_check[1] == "2D info given in parameters but structure is 3D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}".format(True, False, ['INTERFACE', '<NRBASIS>'])


    def test_check_2Dinput_consistency_4(self):
        # case 3: 2D structure but interface parameter set to False
        from aiida_kkr.tools.common_workfunctions import check_2Dinput_consistency
        from aiida.plugins import DataFactory
        StructureData = DataFactory('structure')
        Dict = DataFactory('dict')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        s.set_pbc((True, True, False))
        p = Dict(dict={'INTERFACE':False, '<NRBASIS>':1, '<RBLEFT>':[0,0,0], '<RBRIGHT>':[0,0,0], 'ZPERIODL':[0,0,0], 'ZPERIODR':[0,0,0], '<NLBASIS>':1})
        input_check = check_2Dinput_consistency(s, p)
        assert not input_check[0]
        assert input_check[1] == "'INTERFACE' parameter set to False but structure is 2D"

    def test_check_2Dinput_consistency_5(self):
        # case 5: 3D structure but 2D params given
        from aiida_kkr.tools.common_workfunctions import check_2Dinput_consistency
        from aiida.plugins import DataFactory
        StructureData = DataFactory('structure')
        Dict = DataFactory('dict')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        s.set_pbc((True, True, True))
        p = Dict(dict={'INTERFACE':True, '<NRBASIS>':1, '<RBLEFT>':[0,0,0], '<RBRIGHT>':[0,0,0], 'ZPERIODL':[0,0,0], 'ZPERIODR':[0,0,0], '<NLBASIS>':1})
        input_check = check_2Dinput_consistency(s, p)
        assert not input_check[0]
        assert list(input_check[1]).sort() == list("3D info given in parameters but structure is 2D\nstructure is 2D? {}\ninput has 2D info? {}\nset keys are: {}".format(False, True, ['ZPERIODL', '<NRBASIS>', '<RBLEFT>', 'INTERFACE', '<NLBASIS>', 'ZPERIODR', '<RBRIGHT>'])).sort()


    def test_update_params_wf(self):
        from aiida_kkr.tools.common_workfunctions import update_params_wf
        from masci_tools.io.kkr_params import kkrparams
        from aiida.plugins import DataFactory
        Dict = DataFactory('dict')

        k = kkrparams(LMAX=2)
        node1 = Dict(dict=k.values)
        node2 = Dict(dict={'nodename': 'my_changed_name', 'nodedesc': 'My description text', 'EMIN': -1, 'RMAX': 10.})

        unode = update_params_wf(node1, node1)
        assert unode.get_dict() == node1.get_dict()

        unode = update_params_wf(node1, node2)

        d0 = node1.get_dict()
        for i in list(d0.keys()):
            if d0[i] is None:
                d0.pop(i)

        d1 = unode.get_dict()
        for i in list(d1.keys()):
            if d1[i] is None:
                d1.pop(i)

        l_identical, l_diff = [], []
        for i in list(d0.keys()):
            if i in list(d1.keys()):
                l_identical.append([i, d0[i], d1[i]])
            else:
                l_diff.append([0, i, d0[i]])
        for i in list(d1.keys()):
            if i not in list(d0.keys()):
                l_diff.append([1, i, d1[i]])

        assert l_identical ==  [[u'LMAX', 2, 2]]
        assert l_diff.sort() == [[1, u'RMAX', 10.0], [1, u'EMIN', -1.0]].sort()
        return node1, node2, unode


    def test_neworder_potential_wf(self):
        from numpy import loadtxt
        from aiida.orm import load_node
        from aiida.plugins import DataFactory
        from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
        from aiida.orm.importexport import import_data
        Dict = DataFactory('dict')
        import_data('files/db_dump_kkrflex_create.tar.gz')
        GF_host_calc = load_node('de9b5093-25e7-407e-939e-9282c4431343').outputs
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.retrieved.open('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = Dict(dict=settings_dict)
        startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.remote_folder)
        assert startpot_imp_sfd.get_object_content(startpot_imp_sfd.filename)[::1000] == u'C12807143D556463084.6+55 7D117 9D-87 0+25\n20.70351.75\n0521259.2+491.0-462. 02621D74112D03547T00 4D02116D502 6D39\n96.20261.50941.4944.7+30 98-29 .5-3625D07193.58104D0773D27252285417D341 9.506544D548447094.9+38 91063 54-08 6D28277.60909.98111'


    def test_vca_check(self):
        from aiida_kkr.tools.common_workfunctions import vca_check
        pass


    def test_kick_out_corestates_wf(self):
        from aiida_kkr.tools.common_workfunctions import kick_out_corestates_wf
        pass

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
if __name__=='__main__':
    from aiida import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv()
    from aiida.plugins import DataFactory
    StructureData = DataFactory('structure')
    Dict = DataFactory('dict')

    t = Test_common_workfunctions()

    t1 = t.test_update_params_wf()
    #t2 = t.test_prepare_VCA_structure_wf()
    #t3 = t.test_prepare_2Dcalc_wf()
    #t4 = t.test_test_and_get_codenode()
    #t5 = t.test_get_inputs_kkr()
    #t7 = t.test_get_parent_paranode()
#"""
