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
class Test_common_workfunctions_rmq(object):
    """
    Tests for the common workfunctions from tools.common_workfunctions,
    i.e. functions commonly used in this plugin that depend on aiida stuff to work
    these tests use rabbitMQ
    """

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
        from aiida.tools.importexport import import_data
        Dict = DataFactory('dict')
        import_data('files/db_dump_kkrflex_create.tar.gz')
        GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6').outputs
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.retrieved.open('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = Dict(dict=settings_dict)
        startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.remote_folder)
        assert startpot_imp_sfd.get_object_content(startpot_imp_sfd.filename)[::1000] == u'C12807143D556463084.6+55 7D117 9D-87 0+25\n20.70351.75\n0521259.2+491.0-462. 02621D74112D03547T00 4D02116D502 6D39\n96.20261.50941.4944.7+30 98-29 .5-3625D07193.58104D0773D27252285417D341 9.506544D548447094.9+38 91063 54-08 6D28277.60909.98111'

