# -*- coding: utf-8 -*-
"""
@author: ruess
"""
from __future__ import print_function

from builtins import range
from builtins import object
import pytest

@pytest.mark.usefixtures("aiida_env")
class Test_common_workfunctions(object):
    """
    Tests for the common workfunctions from tools.common_workfunctions,
    i.e. functions commonly used in this plugin that depend on aiida stuff to work
    """
    
    def test_generate_inputcard_from_structure(self):
        from aiida_kkr.tools.common_workfunctions import generate_inputcard_from_structure
        from aiida.orm import DataFactory
        StructureData = DataFactory('structure')
        ParameterData = DataFactory('parameter')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        p = ParameterData(dict={'LMAX':2, 'NSPIN':2, 'RMAX':10, 'GMAX':100})
        generate_inputcard_from_structure(p, s, 'inputcard')
        txt = open('inputcard', 'r').readlines()
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
    
    
    def test_check_2Dinput_consistency(self):
        from aiida_kkr.tools.common_workfunctions import check_2Dinput_consistency
        from aiida.orm import DataFactory
        StructureData = DataFactory('structure')
        ParameterData = DataFactory('parameter')
        s = StructureData(cell=[[0.5, 0.5, 0], [1,0,0], [0,0,1]])
        s.append_atom(position=[0,0,0], symbols='Fe')
        p = ParameterData(dict={'INTERFACE':True})
        input_check = check_2Dinput_consistency(s, p)
        assert input_check[0]
        assert input_check[1] == '2D consistency check complete'
        
        
    def test_update_params_wf(self):
        from aiida_kkr.tools.common_workfunctions import update_params_wf
        from masci_tools.io.kkr_params import kkrparams
        from aiida.orm import DataFactory
        ParameterData = DataFactory('parameter')
        
        k = kkrparams(LMAX=2)
        node1 = ParameterData(dict=k.values)
        node2 = ParameterData(dict={'nodename': 'my_changed_name', 'nodedesc': 'My description text', 'EMIN': -1, 'RMAX': 10.})
        
        unode = update_params_wf(node1, node1)
        assert unode == node1
        
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
        assert l_diff == [[1, u'RMAX', 10.0], [1, u'EMIN', -1.0]]
        return node1, node2, unode


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
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')
    ParameterData = DataFactory('parameter')
    
    t = Test_common_workfunctions()
    
    t1 = t.test_update_params_wf()
    #t2 = t.test_prepare_VCA_structure_wf()
    #t3 = t.test_prepare_2Dcalc_wf()
    #t4 = t.test_test_and_get_codenode()
    #t5 = t.test_get_inputs_kkr()
    #t6 = t.test_get_inputs_voronoi()
    #t7 = t.test_get_parent_paranode()
#"""
    
